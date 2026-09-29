CRC organoid DESeq2 analysis
This repository contains the DESeq2 analysis script and gene-level result tables.
Analysis
Bulk RNA-seq data from patient-derived colorectal cancer organoids were analyzed across five culture conditions and three patients, with one sample for each patient–condition combination.
Differential expression was calculated from raw integer counts using the paired design:
design = ~ patient + condition
Complete medium was used as the reference condition. P values were adjusted using the Benjamini–Hochberg method, with padj < 0.05 used as the statistical significance threshold.
Log2 fold changes were additionally shrunk using apeglm with apeMethod = "nbinomCR".
For downstream interpretation:
•	pvalue and padj were taken from the *_raw.tsv tables;
•	effect sizes were taken from log2FoldChange in the *_lfcShrink_apeglm.tsv tables.
The result tables use Ensembl gene identifiers.
Repository contents
•	scripts/run_deseq2.R - DESeq2 analysis script.
•	results/deseq2/*_raw.tsv - unshrunken DESeq2 results containing Wald statistics, P values and adjusted P values.
•	results/deseq2/*_lfcShrink_apeglm.tsv - results containing apeglm-shrunken log2 fold changes.
Software environment
The supplied results were generated using:
•	R 4.5.1
•	DESeq2 1.48.2
•	apeglm 1.30.0
Small numerical differences in apeglm-shrunken estimates may occur when using different package or runtime versions.

