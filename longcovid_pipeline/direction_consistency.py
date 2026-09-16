"""
direction_consistency.py

Checks whether genes that are statistically significant in MULTIPLE datasets
actually change in the SAME direction (up or down) in each -- a check that is
easy to skip but important not to.

Why this exists
----------------
Finding a gene significant in two independent datasets is often reported as
"replication" -- but significance alone does not confirm the same biological
effect. A gene could be significantly UP in one cohort and significantly DOWN
in another (e.g. due to batch effects, different comparator groups, or
different disease subtypes), which is a very different and much weaker
finding than true concordant replication. This script makes that distinction
explicit and machine-checkable rather than relying on manual inspection.

Usage
-----
    python direction_consistency.py \
        --inputs dataset1.csv dataset2.csv \
        --labels DatasetA DatasetB \
        --gene-col Gene --pvalue-col pvalue --lfc-col log2FC \
        --alpha 0.05 \
        --output overlap_direction_check.csv
"""
import argparse
import pandas as pd
from functools import reduce


def load_sig(path, gene_col, pvalue_col, lfc_col, alpha, sep=',', label=None):
    df = pd.read_csv(path, sep=sep)
    df = df[[gene_col, pvalue_col, lfc_col]].dropna()
    df = df.groupby(gene_col, as_index=False).first()
    sig = df[df[pvalue_col] < alpha].copy()
    sig = sig.rename(columns={lfc_col: f'log2FC__{label}', pvalue_col: f'padj__{label}'})
    return sig[[gene_col, f'log2FC__{label}', f'padj__{label}']].rename(columns={gene_col: 'Gene'})


def check_overlap(inputs, labels, gene_col, pvalue_col, lfc_col, alpha, seps=None):
    seps = seps or [','] * len(inputs)
    sig_frames = [load_sig(p, gene_col, pvalue_col, lfc_col, alpha, sep=s, label=l)
                  for p, s, l in zip(inputs, seps, labels)]
    merged = reduce(lambda l, r: pd.merge(l, r, on='Gene', how='inner'), sig_frames)

    lfc_cols = [c for c in merged.columns if c.startswith('log2FC__')]
    if len(lfc_cols) >= 2:
        signs = merged[lfc_cols].apply(lambda row: (row > 0).nunique() == 1, axis=1)
        merged['direction_consistent'] = signs
    return merged.sort_values('direction_consistent', ascending=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--inputs', nargs='+', required=True, help='Paths to DE result files (2 or more).')
    ap.add_argument('--labels', nargs='+', required=True, help='Short labels for each input, same order (e.g. dataset names).')
    ap.add_argument('--gene-col', default='Gene')
    ap.add_argument('--pvalue-col', default='padj', help="Column to threshold on (default: 'padj').")
    ap.add_argument('--lfc-col', default='log2FC', help="Log2 fold-change column name (default: 'log2FC').")
    ap.add_argument('--alpha', type=float, default=0.05, help='Significance threshold (default: 0.05).')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    if len(args.inputs) != len(args.labels):
        raise SystemExit("--inputs and --labels must have the same number of entries.")

    result = check_overlap(args.inputs, args.labels, args.gene_col, args.pvalue_col,
                            args.lfc_col, args.alpha)
    result.to_csv(args.output, index=False)

    n_total = len(result)
    n_consistent = result['direction_consistent'].sum() if 'direction_consistent' in result.columns else 'n/a'
    print(f"{n_total} genes significant (padj < {args.alpha}) in ALL {len(args.inputs)} provided datasets.")
    print(f"Of these, {n_consistent} changed in the SAME direction in every dataset.")
    print(f"\nFull results written to {args.output}")
    if n_total > 0:
        print("\nDirection-consistent genes (most reliable candidates):")
        cons = result[result['direction_consistent']] if 'direction_consistent' in result.columns else result
        print(cons['Gene'].tolist())


if __name__ == '__main__':
    main()
