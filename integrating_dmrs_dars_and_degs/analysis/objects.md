Objects used in integrative analyses
================
Peter Hickey
4 July 2017

-   [Summary](#summary)
-   [Assays](#assays)
    -   [WGBS](#wgbs)
    -   [ATAC-seq](#atac-seq)
    -   [RNA-seq](#rna-seq)
-   [Features](#features)
    -   [hg19](#hg19)
    -   [GENCODE v19](#gencode-v19)
    -   [CGIs](#cgis)
    -   [Promoters](#promoters)
    -   [Enhancers](#enhancers)
    -   [chromHMM](#chromhmm)
-   [Gene-promoter-enhancer linked data](#gene-promoter-enhancer-linked-data)
-   [Functions](#functions)

Summary
=======

Prepare objects for use in all integrative analyses. These objects are saved in [`../objects/assays-and-features.rda`](../objects/assays-and-features.rda)

Assays
======

WGBS
----

### DMRs

-   `dmrs_pos`: The 13,074 F-stat POS DMRs
    -   `big_dmrs_pos`: The 3,659 / 13,074 F-stat POS DMRs where `abs(maxDiff) > 0.5`
    -   `dmrs_NAvsBA9pos`: The 12,895 / 13,074 F-stat POS DMRs where `NAvsBA9pos == TRUE`
        -   `dmrs_NAvsBA9pos_hypo`: The 3,311 / 12,895 DMRs where NA\_pos is hypomethylated relative to BA9\_pos
        -   `dmrs_NAvsBA9pos_hyper`: The 9,584 / 12,895 DMRs where NA\_pos is hypermethylated relative to BA9\_pos
        -   `big_dmrs_NAvsBA9pos`: The 3,156 / 12,895 POS DMRs where `NAvsBA9pos == TRUE & abs(NAvsBA9pos_meanDiff) > 0.5`
            -   `big_dmrs_NAvsBA9pos_hypo`: The 985 / 3,156 bigDMRs where NA\_pos is hypomethylated relative to BA9\_pos
            -   `big_dmrs_NAvsBA9pos_hyper`: The 2,171 / 3,156 bigDMRs where NA\_pos is hypermethylated relative to BA9\_pos

### DMR-CpGs

-   `cpgs`: The 23,059,530 CpGs used in DMR testing
    -   `dmrs_pos_cpgs`: The 255,537 CpGs within `dmrs_pos`
    -   `dmrs_NAvsBA9pos_cpgs`: The 251,224 CpGs within `dmrs_NAvsBA9pos`
        -   `dmrs_NAvsBA9pos_hypo_cpgs`: The 52,269 CpGs within `dmrs_NAvsBA9pos_hypo`
        -   `dmrs_NAvsBA9pos_hyper_cpgs`: The 198,955 CpGs within `dmrs_NAvsBA9pos_hyper`
    -   `big_dmrs_NAvsBA9pos_cpgs`: The 53,264 CpGs within `big_dmrs_NAvsBA9pos`
        -   `big_dmrs_NAvsBA9pos_hypo_cpgs`: The 12,517 CpGs within `big_dmrs_NAvsBA9pos_hypo`
        -   `big_dmrs_NAvsBA9pos_hyper_cpgs`: The 40,747 CpGs within `big_dmrs_NAvsBA9pos_hyper`

ATAC-seq
--------

-   Open chromatin regions (OCRs) based on MACS2 `narrowPeak` output files
-   Condition-specific: `ocrs_NAcc_pos`, `ocrs_NAcc_neg`, `ocrs_BA9_pos`, `ocrs_BA9_neg`
-   Unions: `ocrs_pos`, `ocrs_neg`, `ocrs_overall`
-   Differentially accessible regions (DARs)
-   `dars_pos`, `dars_neg`, `dars_pos_vs_neg`

RNA-seq
-------

-   `degs`: The 2,952 / 24,161 DEGs (`adj.P.Val < 0.05`) between NA\_pos and BA9\_pos (recall that there are 33,351 genes but we only test 24,161 for DE because the remainder are 'unexpressed' in a large subset of the sample)
    -   `degs_pc`: The 2,402 / 19,823 DEGs (`adj.P.Val < 0.05`) between `NA_pos` and `BA9_pos` that are protein-coding genes
    -   `degs_lnc`: The 550 / 13528 DEGs (`adj.P.Val < 0.05`) between `NA_pos` and `BA9_pos` that are lncRNA genes

Features
========

hg19
----

-   `sl`: The lengths of the autosomes in hg19

GENCODE v19
-----------

-   `union`
    -   `genic`
    -   `promoter`
    -   `five_utr`
    -   `three_utr`
    -   `exonic`
    -   `intronic`
    -   `intergenic`
-   `pc_transcripts`
    -   `genic`
    -   `promoter`
    -   `five_utr`
    -   `three_utr`
    -   `exonic`
    -   `intronic`
    -   `intergenic`
-   `lnc_transcripts`
    -   `genic`
    -   `promoter`
    -   `exonic`
    -   `intronic`
    -   `intergenic`

Plus, all of the above linked to DEGs and non-DEGs.

For reference, the 'union' features have total sizes (bp) of:

``` r
sort(sapply(gencode_features$union, function(x) sum(width(x))), decreasing = TRUE)
#>      genic intergenic   intronic   promoter     exonic  three_utr 
#> 1441985937 1439047349 1384081577  187110443   89835035   37496463 
#>   five_utr 
#>   10360697
```

CGIs
----

-   `cgi`
-   `shores`
-   `shelves`
-   `open_sea`

Promoters
---------

Each promoter enhancer is marked as 'active' if it overlaps an OCR

-   `promoters_by_gene`
    -   `promoters_by_gene_pc`
    -   `promoters_by_gene_lnc`
    -   `cgi_promoters_by_gene`
        -   `cgi_promoters_by_gene_pc`
        -   `cgi_promoters_by_gene_lnc`
    -   `shores_promoters_by_gene`
        -   `shores_promoters_by_gene_pc`
        -   `shores_promoters_by_gene_lnc`

Enhancers
---------

We have three sources of putative enhancers:

1.  `permissive_enhancers`: 43,011 permissive enhancers based on CAGE data from the FANTOM5 consortium (defined in Andersson et al. <http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed>)
2.  `tssa_enhancer_pairs`: 72,019 TSS-associated (`tssa`) enhancer pairs based on CAGE data from the FANTOM5 project (defined in Andersson et al. <http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed>). TSS-links have been made against TSSs of protein-coding genes in GENCODE v19 (see `../../genomic-features/scripts/fantom5-enhancer-tss-associations.R`). **NOTE:** There are 25,683 unique enhancers in this list, however, each enhancer may be linked to multiple TSSs, hence why there are many more TSS-linked enhancer pairs than there are enhancers
3.  `H3K27ac_brain`: 82,041 putative enhancer regions **in brain** based on H3K27ac mark (defined in Vermunt et al.)

We further restrict and overlap these to create:

-   `brain_permissive_enhancers`: The 14,643 / 43,011 permissive enhancers that overlap H3K27ac-based brain enhancer regions
-   `brain_tssa_enhancer_pairs`: The 30,321 / 72,019 TSS-associated enhancer pairs where additionally the enhancer overlaps a H3K27ac-based brain enhancer regions
    -   `deg_brain_tssa_enhancer_pairs`: The 4,581 / 30,321 TSS-associated enhancer pairs where additionally the TSS is of a DEG (NA\_pos vs BA9\_pos)
    -   `non_deg_brain_tssa_enhancers`: The 25,740 / 30,321 TSS-associated enhancer pairs where additionally the TSS is of a non-DEG (NA\_pos vs BA9\_pos)
    -   **NOTE:** `subsetByOverlaps(non_deg_brain_tssa_enhancer_pairs, deg_brain_tssa_enhancer_pairs, invert = TRUE, type = "equal")` returns those TSS-associated enhancer pairs where the enhancer is not linked to any TSS of a DEG (NA\_pos vs BA9\_pos)

Furthermore, each FANTOM5 enhancer is marked as 'active' if it overlaps an OCR in the 'overall' set

For the TSS-associated enhancers, we also contruct 'enhancers-by-gene' objects:

-   `fantom5_enhancers_by_gene`: 11,763 genes linked to at least one FANTOM5 enhancer (median = 3)
    -   `brain_fantom5_enhancers_by_gene`: 8,118 genes linked to at least one enhancer in `brain_fantom5_enhancers` (median = 2)

chromHMM
--------

-   `AH46921`: chromHMM tracks for E068 (Anterior caudate; adjacent to NAcc)
-   `AH46922`: chromHMM tracks for E069 (Cingulate gyrus; BA24 is a subset of this region)
-   `AH46294`: chromHMM track for E071 (Hippocampus middle; HC)
-   `AH46926`: chromHMM tracks for AH46926 (Dorsolateral prefrontal cortex, BA9)

For reference, the chromHMM features have total sizes (Mb) of:

|                            |  AH46921|  AH46922|  AH46924|  AH46926|
|----------------------------|--------:|--------:|--------:|--------:|
| Active TSS                 |     32.0|     25.9|     25.8|     29.5|
| Bivalent Enhancer          |      1.0|      2.0|      2.7|      1.5|
| Bivalent/Poised TSS        |      2.1|      3.4|      3.4|      2.9|
| Enhancers                  |     87.0|     86.1|    101.4|     61.0|
| Flanking Active TSS        |     17.4|     15.5|     19.6|     11.7|
| Flanking Bivalent TSS/Enh  |      0.9|      1.5|      2.8|      1.0|
| Genic enhancers            |      8.9|      8.4|     10.6|      7.7|
| Heterochromatin            |     67.3|     39.0|     60.4|     35.2|
| Quiescent/Low              |   1996.7|   1977.1|   1926.1|   2076.0|
| Repressed PolyComb         |      8.9|     18.5|     16.4|     11.6|
| Strong transcription       |     93.9|     76.7|     90.0|     67.4|
| Transcr. at gene 5' and 3' |      0.8|      0.6|      0.8|      0.7|
| Weak Repressed PolyComb    |    156.6|    241.1|    229.7|    183.1|
| Weak transcription         |    402.9|    381.8|    386.0|    388.6|
| ZNF genes & repeats        |      4.7|      3.5|      5.3|      3.1|

Gene-promoter-enhancer linked data
==================================

We want to understand the relationship between promoter and enhancer DMRs and DARs with the expression of the linked gene. That is, for each gene we tested for differential expression (n = ), we want:

-   the expression estimates in the 11 NA\_pos/BA9\_pos RNA-seq samples
-   an indicator if that gene is differentially expressed
-   the accessibility estimates of the ATAC-seq regions in the 11 NA\_pos/BA9\_pos samples and an indicator if the region is a DAR for regions that:
    -   overlap promoters of that gene
    -   overlap enhancers of that gene
-   the average methylaton level of the DMRs in the 12 NA\_pos/BA9\_pos samples that:
    -   overlap promoters of that gene
    -   overlap enhancers of that gene

**NOTE:** By construction, all methylation measurements are at DMRs

This is stored in `rna_atac_meth`, a data frame with list columns for the nested data

Functions
=========

-   `plotRAM()` for plotting RNA-seq, ATAC-seq, and WGBS data for a gene
