---
output: html_document
title: Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq)
---

## Introduction

This document describes the assay for transposase-accessible chromatin using sequencing (ATAC-seq) genome browser tracks for the project 'Neuronal Brain Region-Specific DNA Methylation and Chromatin Accessibility are Associated with Neuropsychiatric Disease Heritability'.
For an overview of the project, please see the [README](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/docs/README.html).

## Mapping and quality control of ATAC-seq reads

We trimmed reads of their adapter sequences using [**trimadap** (v0.1)](https://github.com/lh3/trimadap/archive/0.1.zip) with the following parameters: `trimadap-mt -3 CTGTCTCTTATACACATCTCCGAGCCCACGAGA ${READ1}; trimadap-mt -3 CTGTCTCTTATACACATCTGACGCTGCCGACGA ${READ2}`

We then aligned these trimmed reads to the hg19 build of the human genome (including autosomes, sex chromosomes, mitochondrial sequence, unplaced sequence, and unlocalized sequence) using [**Bowtie2** (v2.2.5)](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) with alignment parameters: `bowtie2 -X 2000 --local --dovetail`. 
Potential PCR duplicate reads were marked using MarkDuplicates from the [Picard library (v2.2.1)](http://broadinstitute.github.io/picard/).

## Identifying differentially accessible regions (DARs)

Peaks were called in each condition (NAcc_pos, NAcc_neg, BA9_pos, and BA9_neg) using [MACS (v2.1.0)](https://github.com/taoliu/MACS). Specifically, data from each condition were combined into a metasample formed by using all non-duplicate-marked reads with a mapping quality > 30 and then processed using: `macs2 callpeaks --nomodel --nolambda --call-summits -t ${BAMS[@]}`. 
We took the 'narrowPeaks' produced by MACS and filtered out those regions overlapping the [ENCODE mappability consensus blacklist regions](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/) and the [blacklist for ATAC-seq created by Buenrostro et al (access requires signing in with a Google Account)](https://sites.google.com/site/atacseqpublic/atac-seq-analysis-methods/mitochondrialblacklists-1). 
We took this filtered list as our condition-specific sets of open chromatin regions (OCRs).

To perform the differential analysis, we first took the union of condition-specific OCRs on the autosomes to construct an 'overall' set of OCRs. 
This 'overall' set of OCRs contained 853,053 regions (630 Mb). 
For each sample, we then counted the number of fragments (fragment = start of read1 to end of read2) overlapping each of the 'overall' OCRs using the `summarizeOverlaps()` function in the [**GenomicAlignments** R/Bioconductor package (v1.10.0)](https://bioconductor.org/packages/GenomicAlignments/). 
Specifically, we only counted those fragments where both reads had a mapping-quality score > 30, reads not marked as potential PCR duplicates, and those where any part of the fragment overlapped exactly one peak.

We then analyzed these data using the [**voom** method](https://www.bioconductor.org/packages/limma/), originally designed for differential expression analysis of RNA-seq data17. 
Briefly, the read counts were transformed to counts per million (cpm) and only those 283,812 / 853,053 peaks with 
at least 1 cpm for at least 5 samples (the size of the smallest group of samples) were retained. 
We normalized these counts using TMM, then used [**edgeR** (v3.16.5)](https://www.bioconductor.org/packages/edgeR/) and [**limma** (v3.30.7)](https://www.bioconductor.org/packages/limma/) to transform these counts to log2-cpm, estimate the mean-variance relationship, and compute appropriate observation-level weights ready for linear modelling.

In our design matrix, we blocked on donor (donor1, ..., donor6) and included a term for each group (e.g., BA9_neg for NeuN- cells from BA9, BA9_pos for NeuN+ cells from BA9, etc.). 
We ran surrogate variable analysis21 using the [**sva** (v3.22.0) R/Bioconductor package](https://www.bioconductor.org/packages/sva/) and identified 4 surrogate variables, one of which correlated with the date on which these samples were flow-sorted. 
We ultimately decided to include all 4 surrogate variables in the linear model. 
Using the empirical Bayes shrinkage method implemented in [**limma**](https://www.bioconductor.org/packages/limma/), we tested for differential accessibility of peaks in three comparisons: 

1. NAcc vs. BA9 in NeuN+ cells
2. NAcc vs. BA9 in NeuN- cells
3. NeuN+ cells vs. NeuN- cells

For an ATAC-seq peak to be called a differentially accessible region (DAR), it had to have a Benjamini-Hochberg adjusted P-value < 0.05.

## Credits

Questions should be directed to [Peter Hickey](mailto:peter.hickey@gmail.com).
