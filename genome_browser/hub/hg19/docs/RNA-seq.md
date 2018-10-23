---
output: html_document
title: RNA-seq
---

## Introduction

This document describes the RNA-seq genome browser tracks for the project 'Neuronal Brain Region-Specific DNA Methylation and Chromatin Accessibility are Associated with Neuropsychiatric Disease Heritability'.
For an overview of the project, please see the [README](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/docs/README.html).


## Mapping and quality control of RNA-seq reads

We trimmed the first 3 bp of read1, which were derived from template switching oligos and not the cDNA of interest, using [**seqtk** (v1.2-r94)](https://github.com/lh3/seqtk) with the following parameters:
`seqtk trimfq -b 3 ${READ1}`. 
We then quasi-mapped these trimmed reads to a FASTA file of protein-coding and lncRNA genes from [GENCODE v19](http://www.gencodegenes.org/releases/19.html) and performed transcript-level quantification using [**Salmon** (v0.7.2)](https://github.com/COMBINE-lab/salmon).

## Identifying differentially expressed genes (DEGs)

We used [**tximport** (v1.2.0)](https://www.bioconductor.org/packages/tximport/) to compute normalized gene-level counts from the transcript-level abundance estimates (scaling these using the average transcript length over samples and the library size). 
Only autosomal genes with at least 1 cpm in at least 4 libraries (the size of the smallest group of samples) were retained for downstream analysis (24,161 / 33,351 genes). 
We normalized these counts using TMM18 then used [**edgeR** (v3.16.5)](https://www.bioconductor.org/packages/edgeR/) and [**limma** (v3.30.7)](https://www.bioconductor.org/packages/limma/) to transform these counts to log2-cpm, estimate the mean-variance relationship, and compute appropriate observation-level weights ready for linear modelling.

In our design matrix, we blocked on donor (donor1, ..., donor6) and included a term for each group (e.g., BA9_neg for NeuN- cells from BA9, BA9_pos for NeuN+ cells from BA9, etc.). 
We ran surrogate variable analysis21 using the [**sva** (v3.22.0) R/Bioconductor package](https://www.bioconductor.org/packages/sva/) and identified 5 surrogate variables, some of which correlated with the date on which these samples were flow-sorted. 
We ultimately decided to include all 5 surrogate variables in the linear model. 
Using the empirical Bayes shrinkage method implemented in [**limma**](https://www.bioconductor.org/packages/limma/), we tested for differential expression of genes in three comparisons:

1. NAcc vs. BA9 in NeuN+ cells
2. NAcc vs. BA9 in NeuN- cells
3. NeuN+ cells vs. NeuN- cells

For a gene to be called a differentially expressed gene (DEG), it had to have a Benjamini-Hochberg adjusted P-value < 0.05 with no minimum log2 fold change cutoff.

## Credits

Questions should be directed to [Peter Hickey](mailto:peter.hickey@gmail.com).
