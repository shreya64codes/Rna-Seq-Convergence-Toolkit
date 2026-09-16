"""
longcovid_pipeline

A small, reusable toolkit for harmonized cross-dataset differential
expression analysis, originally developed for a Long COVID transcriptomics
meta-analysis project but applicable to any disease area combining multiple
small, independent public RNA-seq/microarray datasets.

Modules
-------
prepare_for_deseq2   -- reshape a combined counts matrix into Galaxy-DESeq2-ready
                         per-sample files + design table
rank_aggregation     -- combine DE results across datasets via percentile rank
direction_consistency-- check that genes significant in multiple datasets
                         actually change in the same direction
specificity_check    -- test whether candidate genes are specific to a
                         condition or shared with a comparator condition
"""
