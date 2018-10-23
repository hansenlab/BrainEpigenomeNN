---
output: html_document
title: Neuronal brain region-specific DNA methylation and chromatin accessibility are associated with neuropsychiatric disease heritability
---

## Description

These tracks show results from genome-wide analyses of DNA methylation, chromatin accessibility, and gene expression data from 15 human brain donors [[1](https://www.biorxiv.org/content/early/2017/03/24/120386)].

We performed:

- **W**hole **G**enome **B**isulfite **S**equencing (WGBS, n = 72)
- **A**ssay for **T**ransposase-**A**ccessible **C**hromatin using sequencing 
(ATAC-seq, n = 22)
- RNA sequencing (RNA-seq, n = 20) 

Data come from 4 brain regions:

- Dorsolateral prefrontal cortex (BA9)
- Anterior cingulate gyrus (BA24)
- Hippocampus (HC)
- Nucleus accumbens (NAcc) 

and 2 'cell types' and bulk tissue:

- NeuN+ cells isolated using fluourescence activated nuclei sorting (NeuN+)
- NeuN- cells isolated using fluourescence activated nuclei sorting (NeuN-)
- Bulk tissue (Unsorted)

<img src="https://s3.us-east-2.amazonaws.com/brainepigenome/img/individual-level_dataset_matrix.png" alt="Individual-level dataset matrix" style="width: 600px;"/>

<img src="https://s3.us-east-2.amazonaws.com/brainepigenome/img/condition-level_dataset_matrix.png" alt="Condition-level dataset matrix" style="width: 600px;"/>

## Methods

Methods for each assay are described in the following documents:

- [**W**hole **G**enome **B**isulfite **S**equencing (WGBS)](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/docs/WGBS.html)
- [**A**ssay for **T**ransposase-**A**ccessible **C**hromatin using sequencing (ATAC-seq)](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/docs/ATAC-seq.html)
- [RNA sequencing (RNA-seq)](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/docs/RNA-seq.html)

Please refer to [[1](https://www.biorxiv.org/content/early/2017/03/24/120386)] for further details.

## Caveats

Methylation and coverage tracks do not incorporate biological variability. 
Visually inferring regions of differential methylation, differential accessibility, or differential expression from these tracks alone is not recommended.
Instead, information across all samples should be used to identify these regions and analysed with statistical methods that account for biological variability (see 'Data Access', below). 

## Data Access

Raw and processed data generated during this study will be made available on the Gene Expression Omnibus ([GSE96615](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96615)) upon publication. 
We recommend that these data are used for any downstream statistical analyses.

The raw data for these tracks can be explored interactively using the [Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables) or [Data Integrator](http://genome.ucsc.edu/cgi-bin/hgIntegrator). 
The tables can also be downloaded from our downloads server for local processing; the exact filenames can be found in the [track configuration file](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/trackDb.txt) under the `bigDataUrl` field. To download a track, append the `bigDataUrl` to the base URL of [https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/)
For example, the track 'Average mCG (small smooth) in BA9 (NeuN+) samples' has `bigDataUrl` of `BA9_pos.small_smooth.mCG.bw` and can be downloaded from [https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/BA9_pos.small_smooth.mCG.bw](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/BA9_pos.small_smooth.mCG.bw).

## Credits

Questions should be directed to [Peter Hickey](mailto:peter.hickey@gmail.com).

## References

1. [Rizzardi, L. et al. Neuronal brain region-specific DNA methylation and chromatin accessibility are associated with neuropsychiatric disease heritability. bioRxiv 120386 (2017). doi:10.1101/120386](https://www.biorxiv.org/content/early/2017/03/24/120386)
