# Long COVID Cross-Dataset Convergence Pipeline

A small, reusable Python toolkit for **harmonized differential expression
meta-analysis across multiple small, independent public datasets**.

This pipeline was built during a project reanalyzing four independent public
Long COVID transcriptomic datasets from NCBI GEO, but every script here is
disease-agnostic — it works for any situation where you have several small
RNA-seq/microarray cohorts and want to find candidate genes that hold up
across studies, rather than trusting any single underpowered cohort alone.

## Why this exists

Individual small cohorts (common in rare diseases, post-viral conditions,
and emerging conditions like Long COVID) are usually **statistically
underpowered**: after correcting for testing thousands of genes at once,
very few genes survive strict significance in any one study, even when real
biological signal is present. Meanwhile, naively combining raw data from
different studies is unsafe, because:

- studies use different **case definitions**, **comparator groups**, and
  **sample types** (e.g. PBMC vs. whole blood)
- **repeated/longitudinal samples** from the same patient can masquerade as
  independent replication and inflate significance
- genes can be **significant in multiple datasets but change in opposite
  directions** — apparent "replication" that isn't real
- a **general infection/inflammation response** can be mistaken for a
  condition-specific signal unless checked against a comparator condition

This pipeline provides small, focused tools to handle each of these issues
explicitly rather than silently.

## Pipeline overview

```
Raw counts matrix + metadata (per dataset)
        │
        ▼
[1] prepare_for_deseq2.py   ── reshape into Galaxy-DESeq2-ready format
        │                       (also handles patient deduplication)
        ▼
   Run DESeq2 externally (usegalaxy.org, R/Bioconductor, or pyDESeq2)
        │
        ▼
[2] rank_aggregation.py     ── combine DE results across datasets via
        │                       percentile rank (works even if no single
        │                       dataset reaches strict significance)
        ▼
[3] direction_consistency.py── for genes significant in multiple datasets,
        │                       verify they change in the SAME direction
        ▼
[4] specificity_check.py    ── test candidates against a comparator
                                condition to rule out general/non-specific
                                signal
```

Steps 2-4 are independent and can be run in any order/combination depending
on what you need.

## Why DESeq2 runs externally

This pipeline deliberately does **not** attempt to reimplement DESeq2's
statistics in pure Python. DESeq2's negative-binomial model is the field
standard for RNA-seq count data and is what most published studies use — a
simpler substitute (e.g. a t-test) was tried during this project's
development and produced results that failed to replicate the original
authors' own published findings on the same data. Rather than risk the same
mistake for other users, this pipeline treats DESeq2 as an external step:

- **Recommended, no install required:** [usegalaxy.org](https://usegalaxy.org) —
  free, runs real DESeq2, works entirely in-browser.
  `prepare_for_deseq2.py` outputs are built specifically for Galaxy's DESeq2
  tool (see its docstring for the exact tool settings to use).
- **Alternative:** R/Bioconductor's `DESeq2` package, or the Python port
  [`pydeseq2`](https://pydeseq2.readthedocs.io) if you have an environment
  with internet access to install it.

Once you have DESeq2 output (a table with gene ID, log2FoldChange, pvalue,
padj), everything downstream (`rank_aggregation.py`,
`direction_consistency.py`, `specificity_check.py`) is pure Python with no
external dependencies beyond pandas/numpy/scipy.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### 1. Prepare a dataset for DESeq2

```bash
python longcovid_pipeline/prepare_for_deseq2.py \
    --counts my_dataset_counts.csv \
    --design my_dataset_design.csv \
    --sample-col sample_id \
    --condition-col disease_status \
    --dedup-patient-col patient_id \
    --dedup-sort-col collection_date \
    --outdir prepped/my_dataset
```

This writes `prepped/my_dataset/per_sample/*.tsv` and
`prepped/my_dataset/design_table.tsv`, ready to upload to Galaxy. See the
script's docstring for the exact Galaxy DESeq2 tool settings to use with
these files.

`--dedup-patient-col` / `--dedup-sort-col` are optional — use them whenever
your samples include repeated measurements from the same individuals, to
avoid one of the most common (and easy to miss) sources of falsely inflated
significance.

### 2. Aggregate ranks across datasets

```bash
python longcovid_pipeline/rank_aggregation.py \
    --inputs dataset1_deseq2.csv dataset2_deseq2.csv dataset3_deseq2.csv \
    --gene-col Gene --pvalue-col pvalue \
    --min-datasets 2 \
    --output ranked_all_genes.csv \
    --top-output ranked_top_candidates.csv
```

### 3. Check direction consistency

```bash
python longcovid_pipeline/direction_consistency.py \
    --inputs dataset1_deseq2.csv dataset2_deseq2.csv \
    --labels Dataset1 Dataset2 \
    --gene-col Gene --pvalue-col padj --lfc-col log2FC \
    --output direction_check.csv
```

### 4. Check specificity against a comparator condition

```bash
python longcovid_pipeline/specificity_check.py \
    --candidates GENE1 GENE2 GENE3 \
    --comparator-results acute_condition_deseq2.csv \
    --gene-col Gene --pvalue-col padj \
    --output specificity_results.csv
```

Or from a file of candidate genes:

```bash
python longcovid_pipeline/specificity_check.py \
    --candidates-file my_candidates.txt \
    --comparator-results acute_condition_deseq2.csv \
    --output specificity_results.csv
```

## Worked example

`example_data/` contains small synthetic datasets (2 datasets, 100 genes,
~10 samples each, with a known planted signal in GENE1-GENE3) that
demonstrate the full pipeline end-to-end. Run the commands above against
these files to see expected output before applying the pipeline to your own
data.

## Background: the project this came from

This pipeline was developed while reanalyzing four independent public Long
COVID datasets (GEO accessions GSE224615, GSE270045, GSE275334, GSE226260)
to test whether independent cohorts converge on shared candidate biomarkers.
Key lessons from that project, now encoded as pipeline safeguards:

- An initial t-test-based analysis of one dataset produced results that did
  not replicate the original authors' own DESeq2-based findings — resolved
  by switching to true DESeq2 (`prepare_for_deseq2.py` exists because of
  this).
- One dataset's naive analysis showed an implausible 3,205 "significant"
  genes, traced to repeated longitudinal patient samples being treated as
  independent (`--dedup-patient-col` exists because of this).
- Of 19 genes significant in two independent datasets, only 3 changed in
  the same direction — the rest were a statistical coincidence, not
  replication (`direction_consistency.py` exists because of this).
- Several candidate genes lost their apparent Long-COVID specificity once
  checked against an acute-COVID comparator dataset
  (`specificity_check.py` exists because of this).

Full project write-up, literature review, and results are documented
separately in the accompanying project report.

## License

MIT — reuse freely, with attribution appreciated.
