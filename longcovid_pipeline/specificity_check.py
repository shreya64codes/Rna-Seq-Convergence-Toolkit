"""
specificity_check.py

Screens a list of candidate genes against a comparator/control dataset's
differential-expression results, to test whether a candidate is specific to
the condition of interest or is instead part of a more general biological
response shared with a related condition.

Example use case (this pipeline's origin)
-------------------------------------------
Candidate genes were found differentially expressed between Long COVID
patients and recovered/healthy controls. But some of those genes might just
reflect "having been infected with COVID at all" rather than anything
specific to *failing to recover*. To test this, the same candidate genes were
checked against an ACUTE COVID-19 dataset (active infection vs. non-COVID
controls). Genes that were ALSO significant during acute infection were
reclassified as general infection/inflammation markers rather than
condition-specific candidates.

This script generalizes that pattern to any two-condition comparison.

Usage
-----
    python specificity_check.py \
        --candidates ALOX15B NFASC INTU SRCAP DDX43 \
        --comparator-results acute_condition_deseq2.csv \
        --gene-col Gene --pvalue-col padj --alpha 0.05 \
        --output specificity_results.csv

--candidates can also be a path to a text file with one gene per line:
    python specificity_check.py --candidates-file candidates.txt ...
"""
import argparse
import pandas as pd


def run_specificity_check(candidates, comparator_path, gene_col, pvalue_col, alpha, sep=','):
    comp = pd.read_csv(comparator_path, sep=sep)
    if gene_col not in comp.columns or pvalue_col not in comp.columns:
        raise ValueError(f"Comparator file missing expected columns '{gene_col}'/'{pvalue_col}'. "
                          f"Found: {list(comp.columns)}")

    rows = []
    for gene in candidates:
        match = comp[comp[gene_col] == gene]
        if len(match) == 0:
            rows.append({'Gene': gene, 'found_in_comparator': False,
                         'comparator_pvalue': None, 'specific': None})
            continue
        pval = match.iloc[0][pvalue_col]
        rows.append({
            'Gene': gene,
            'found_in_comparator': True,
            'comparator_pvalue': pval,
            'specific': not (pval < alpha)  # specific = NOT significant in the comparator
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--candidates', nargs='+', help='Candidate gene names, space-separated.')
    group.add_argument('--candidates-file', help='Path to a text file, one gene per line.')
    ap.add_argument('--comparator-results', required=True,
                     help='DE results file for the comparator condition (e.g. acute disease vs. its own controls).')
    ap.add_argument('--gene-col', default='Gene')
    ap.add_argument('--pvalue-col', default='padj')
    ap.add_argument('--alpha', type=float, default=0.05)
    ap.add_argument('--sep', default=',', help="Delimiter for the comparator results file (default: ',').")
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    if args.candidates_file:
        with open(args.candidates_file) as f:
            candidates = [line.strip() for line in f if line.strip()]
    else:
        candidates = args.candidates

    result = run_specificity_check(candidates, args.comparator_results, args.gene_col,
                                    args.pvalue_col, args.alpha, args.sep)
    result.to_csv(args.output, index=False)

    n_specific = result['specific'].sum()
    n_general = (result['specific'] == False).sum()
    n_missing = (~result['found_in_comparator']).sum()

    print(f"Checked {len(candidates)} candidate genes against comparator dataset.")
    print(f"  Specific (not significant in comparator): {n_specific}")
    print(f"  General/shared (also significant in comparator): {n_general}")
    print(f"  Not found in comparator dataset: {n_missing}")
    print(f"\nFull results written to {args.output}")


if __name__ == '__main__':
    main()
