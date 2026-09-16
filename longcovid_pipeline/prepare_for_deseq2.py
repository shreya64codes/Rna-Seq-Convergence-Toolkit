"""
prepare_for_deseq2.py

Converts a combined gene-expression counts matrix (genes x samples) into the
per-sample file format required by Galaxy's DESeq2 tool ("Select contrasts
from sample sheet columns" mode), plus a matching tab-separated design table.

Why this exists
----------------
Galaxy's official DESeq2 wrapper expects a *collection* of individual
per-sample count files (two columns: gene, count; no header), rather than one
merged matrix. Most public datasets (e.g. GEO) are distributed as a single
merged matrix. This script bridges that gap.

Usage
-----
    python prepare_for_deseq2.py \
        --counts counts_matrix.csv \
        --design design.csv \
        --sample-col sample \
        --condition-col condition \
        --outdir output/

Input requirements
-------------------
counts_matrix.csv : genes as rows, samples as columns, first column = gene ID.
                     Values should be raw/estimated integer counts, NOT
                     normalized (TPM/FPKM/CPM) values -- DESeq2 requires raw
                     counts, since it models count-based noise internally.

design.csv         : one row per sample, at minimum two columns:
                        - a sample-identifier column matching the counts
                          matrix's column headers
                        - a condition/group column (e.g. "LC" / "HC")

Output
------
outdir/per_sample/<sample_name>.tsv   -- one two-column file per sample
outdir/design_table.tsv               -- tab-separated design table for Galaxy
                                          (Galaxy's "tabular" format expects
                                          tabs, not commas -- a common gotcha)

Notes on gotchas encountered in practice
-----------------------------------------
- Comma-separated design tables are silently misread by Galaxy as a single
  column. Always emit tab-separated design tables (this script does this
  automatically).
- Repeated/longitudinal samples from the same patient must be deduplicated
  (e.g. keep only the first timepoint) BEFORE running DESeq2, or patient
  identity gets treated as independent replication and inflates significance.
  Use `dedup_by_patient()` below if your design table has a patient ID column.
"""
import argparse
import os
import sys
import pandas as pd


def dedup_by_patient(design: pd.DataFrame, patient_col: str, sort_col: str = None) -> pd.DataFrame:
    """
    Keep only one row per unique patient (default: first occurrence).

    Use this BEFORE prepare_files() whenever samples include repeated
    measurements from the same individual over time. Treating repeated
    samples as independent observations is a classic source of falsely
    inflated significance in differential expression testing.

    Parameters
    ----------
    design : DataFrame with at least a patient identifier column.
    patient_col : column name identifying unique individuals.
    sort_col : optional column to sort by before deduplicating (e.g. a date
               column), so "first occurrence" means "earliest timepoint".
    """
    d = design.copy()
    if sort_col and sort_col in d.columns:
        d[sort_col] = pd.to_datetime(d[sort_col], errors='coerce')
        d = d.sort_values(sort_col)
    return d.drop_duplicates(subset=patient_col, keep='first')


def prepare_files(counts_path: str, design_path: str, sample_col: str,
                   condition_col: str, outdir: str, sep_counts: str = ',',
                   sep_design: str = ','):
    """
    Main entry point: reads a counts matrix + design table, writes per-sample
    count files and a Galaxy-compatible tab-separated design table.
    """
    counts = pd.read_csv(counts_path, sep=sep_counts, index_col=0, low_memory=False)
    design = pd.read_csv(design_path, sep=sep_design)

    if sample_col not in design.columns:
        sys.exit(f"ERROR: sample column '{sample_col}' not found in design table. "
                  f"Available columns: {list(design.columns)}")
    if condition_col not in design.columns:
        sys.exit(f"ERROR: condition column '{condition_col}' not found in design table. "
                  f"Available columns: {list(design.columns)}")

    samples = design[sample_col].tolist()
    missing = [s for s in samples if s not in counts.columns]
    if missing:
        print(f"WARNING: {len(missing)} samples in design table not found in counts "
              f"matrix columns and will be skipped: {missing[:5]}{'...' if len(missing) > 5 else ''}")
        samples = [s for s in samples if s in counts.columns]
        design = design[design[sample_col].isin(samples)]

    persample_dir = os.path.join(outdir, 'per_sample')
    os.makedirs(persample_dir, exist_ok=True)

    counts_subset = counts[samples].apply(pd.to_numeric, errors='coerce').round(0).dropna(how='any').astype(int)

    for sample in samples:
        out = counts_subset[[sample]].reset_index()
        out.columns = ['Gene', 'Count']
        out.to_csv(os.path.join(persample_dir, f'{sample}.tsv'), sep='\t', index=False, header=False)

    design_out = design[[sample_col, condition_col]].rename(
        columns={sample_col: 'sample', condition_col: 'condition'})
    design_path_out = os.path.join(outdir, 'design_table.tsv')
    design_out.to_csv(design_path_out, sep='\t', index=False)

    print(f"Wrote {len(samples)} per-sample files to {persample_dir}")
    print(f"Wrote design table to {design_path_out}")
    print(f"Group sizes: {design_out['condition'].value_counts().to_dict()}")
    print()
    print("Next steps (Galaxy/usegalaxy.org):")
    print("  1. Upload all files in per_sample/ plus design_table.tsv")
    print("  2. Select all per-sample files -> build a Dataset List collection")
    print("  3. Run the DESeq2 tool with 'how' = 'Select contrasts from sample sheet columns'")
    print("  4. Set Count file(s) collection = your new collection")
    print("  5. Set Sample sheet file = design_table.tsv")
    print("  6. Set Factor levels column(s) = condition; set reference/target levels")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--counts', required=True, help='Path to combined counts matrix (genes x samples).')
    ap.add_argument('--design', required=True, help='Path to design/metadata table (one row per sample).')
    ap.add_argument('--sample-col', required=True, help='Column in design table matching counts matrix column headers.')
    ap.add_argument('--condition-col', required=True, help='Column in design table giving the group/condition label.')
    ap.add_argument('--outdir', required=True, help='Output directory.')
    ap.add_argument('--sep-counts', default=',', help="Delimiter for counts file (default: ',').")
    ap.add_argument('--sep-design', default=',', help="Delimiter for design file (default: ',').")
    ap.add_argument('--dedup-patient-col', default=None,
                     help='Optional: column identifying unique patients, to deduplicate repeated/longitudinal samples before analysis.')
    ap.add_argument('--dedup-sort-col', default=None,
                     help='Optional: column to sort by (e.g. a date) so deduplication keeps the earliest timepoint.')
    args = ap.parse_args()

    if args.dedup_patient_col:
        design = pd.read_csv(args.design, sep=args.sep_design)
        design = dedup_by_patient(design, args.dedup_patient_col, args.dedup_sort_col)
        tmp_design_path = os.path.join(os.path.dirname(args.design) or '.', '_deduped_design.csv')
        design.to_csv(tmp_design_path, index=False)
        print(f"Deduplicated design table ({len(design)} unique patients) written to {tmp_design_path}")
        args.design = tmp_design_path
        args.sep_design = ','

    os.makedirs(args.outdir, exist_ok=True)
    prepare_files(args.counts, args.design, args.sample_col, args.condition_col,
                  args.outdir, args.sep_counts, args.sep_design)


if __name__ == '__main__':
    main()
