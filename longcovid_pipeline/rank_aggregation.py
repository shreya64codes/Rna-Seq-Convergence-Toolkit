"""
rank_aggregation.py

Combines DESeq2 (or any differential-expression) result files from multiple
independent datasets into a single ranked candidate gene list, using
percentile-rank averaging rather than requiring strict significance in every
individual dataset.

Why this exists
----------------
Small cohorts are common in rare-disease and post-viral-condition research.
A gene may fail strict significance (padj < 0.05) in every individual small
dataset simply due to low statistical power, while still consistently
ranking near the top of the "most different" gene list across datasets --
which is itself meaningful evidence. This script formalizes that comparison.

Method
------
1. Within each dataset, every gene is ranked by its p-value (rank 1 = most
   significant).
2. Ranks are converted to PERCENTILES (0 to 1), since datasets typically test
   different numbers of genes (a targeted 785-gene panel vs. a 24,000-gene
   whole-transcriptome dataset cannot be compared on raw rank alone).
3. For each gene, percentile ranks are averaged across every dataset in which
   it was measured.
4. Only genes present in at least `--min-datasets` datasets are reported, to
   avoid rewarding a gene simply for appearing in only one study.

Usage
-----
    python rank_aggregation.py \
        --inputs dataset1.csv dataset2.csv dataset3.csv \
        --gene-col Gene --pvalue-col pvalue \
        --min-datasets 2 \
        --output ranked_candidates.csv

Each input file must be a CSV/TSV with at least a gene identifier column and
a p-value column (column names are configurable and must match across files,
or be aligned beforehand).
"""
import argparse
import os
import pandas as pd
from functools import reduce


def load_and_rank(path: str, gene_col: str, pvalue_col: str, sep: str = ',', label: str = None) -> pd.DataFrame:
    """Load one DE results file and compute percentile ranks by p-value."""
    df = pd.read_csv(path, sep=sep)
    if gene_col not in df.columns or pvalue_col not in df.columns:
        raise ValueError(f"{path}: expected columns '{gene_col}' and '{pvalue_col}', "
                          f"found {list(df.columns)}")
    df = df[[gene_col, pvalue_col]].dropna()
    df = df.groupby(gene_col, as_index=False)[pvalue_col].min()  # dedup, keep best p-value
    label = label or os.path.splitext(os.path.basename(path))[0]
    df[f'rank_pct__{label}'] = df[pvalue_col].rank(pct=True)
    return df[[gene_col, f'rank_pct__{label}']].rename(columns={gene_col: 'Gene'})


def aggregate(inputs, gene_col, pvalue_col, seps=None, labels=None):
    """Merge multiple ranked datasets and compute the average percentile rank per gene."""
    seps = seps or [','] * len(inputs)
    labels = labels or [None] * len(inputs)
    ranked = [load_and_rank(p, gene_col, pvalue_col, sep=s, label=l)
              for p, s, l in zip(inputs, seps, labels)]
    merged = reduce(lambda l, r: pd.merge(l, r, on='Gene', how='outer'), ranked)

    rank_cols = [c for c in merged.columns if c.startswith('rank_pct__')]
    merged['n_datasets'] = merged[rank_cols].notna().sum(axis=1)
    merged['avg_rank_pct'] = merged[rank_cols].mean(axis=1)
    return merged.sort_values('avg_rank_pct')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--inputs', nargs='+', required=True, help='Paths to DE result files (one per dataset).')
    ap.add_argument('--gene-col', default='Gene', help="Gene identifier column name (default: 'Gene').")
    ap.add_argument('--pvalue-col', default='pvalue', help="P-value column name (default: 'pvalue').")
    ap.add_argument('--min-datasets', type=int, default=2,
                     help='Minimum number of datasets a gene must appear in to be included in the top-candidates output (default: 2).')
    ap.add_argument('--output', required=True, help='Path to write the full merged ranking (CSV).')
    ap.add_argument('--top-output', default=None,
                     help='Optional: path to write just the top-candidates table (genes meeting --min-datasets).')
    args = ap.parse_args()

    merged = aggregate(args.inputs, args.gene_col, args.pvalue_col)
    merged.to_csv(args.output, index=False)
    print(f"Wrote full merged ranking ({len(merged)} genes) to {args.output}")

    top = merged[merged['n_datasets'] >= args.min_datasets]
    print(f"{len(top)} genes present in >= {args.min_datasets} datasets.")
    print()
    print("Top 15 candidates by average rank percentile:")
    print(top.head(15).to_string(index=False))

    if args.top_output:
        top.to_csv(args.top_output, index=False)
        print(f"\nWrote top-candidates table to {args.top_output}")


if __name__ == '__main__':
    main()
