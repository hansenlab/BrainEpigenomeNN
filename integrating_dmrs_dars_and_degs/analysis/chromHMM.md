Context and enrichment/depetion of DMRs and DAPs with respect to chromHMM features
================
Peter Hickey
9 December 2016

-   [chromHMM context of DMR-CpGs and bigDMR-CpGs](#chromhmm-context-of-dmr-cpgs-and-bigdmr-cpgs)
    -   [DMR-CpGs](#dmr-cpgs)
    -   [bigDMR-CpGs](#bigdmr-cpgs)
    -   [Summary](#summary)
-   [chromHMM context of ATAC peaks, DAPs, and bigDAPs](#chromhmm-context-of-atac-peaks-daps-and-bigdaps)
    -   [ATAC-peaks](#atac-peaks)
    -   [DAPs](#daps)
    -   [bigDAPs](#bigdaps)
    -   [Summary](#summary-1)
-   [chromHMM enrichment/depletion of DMR-CpGs and bigDMR-CpGs](#chromhmm-enrichmentdepletion-of-dmr-cpgs-and-bigdmr-cpgs)
    -   [POS DMRs](#pos-dmrs-3)
    -   [NA\_pos vs. BA9\_pos DMRs](#na_pos-vs.-ba9_pos-dmrs-1)
    -   [Summary](#summary-2)
-   [chromHMM enrichment/depletion of ATAC-seq peaks, DAPs, and bigDAPs](#chromhmm-enrichmentdepletion-of-atac-seq-peaks-daps-and-bigdaps)
    -   [Counting peaks](#counting-peaks)
    -   [Counting bases](#counting-bases)

chromHMM context of DMR-CpGs and bigDMR-CpGs
============================================

Unlike the GENCODE-based features, chromHMM features are cell type-specific. We use 4 different chromHMM tracks:

|     | AH      | E    | Region                       | Info                            |
|-----|:--------|:-----|:-----------------------------|:--------------------------------|
| 2   | AH46921 | E068 | AnteriorCaudate              | Adjacent to Nacc                |
| 3   | AH46922 | E069 | CingulateGyrus               | BA24 is a subset of this region |
| 4   | AH46924 | E071 | HippocampusMiddle            | HC                              |
| 6   | AH46926 | E073 | DorsolateralPrefrontalCortex | BA9                             |

The four tracks are approximately 80% pairwise identical, but there is still a reasonable amount that is unique to each track:

``` r
z <- list("AH46921" = AH46921, "AH46922" = AH46922, "AH46924" = AH46924, 
          "AH46926" = AH46926)

as.data.frame(
  lapply(list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4)), function(i) {
    zz <- GenomicRanges::intersect(z[[i[1]]], z[[i[2]]])
    val <- list(round(sum(as.numeric(sum(width(zz)))) / 
                        sum(as.numeric(sum(width(z[[i[1]]])))) * 100, 0),
                round(sum(as.numeric(sum(width(zz)))) / 
                        sum(as.numeric(sum(width(z[[i[2]]])))) * 100, 0))
    setNames(val, names(z)[i])
  }))
#>   AH46921 AH46922 AH46921.1 AH46924 AH46921.2 AH46926 AH46922.1 AH46924.1
#> 1      81      81        81      81        81      81        82        82
#>   AH46922.2 AH46926.1 AH46924.2 AH46926.2
#> 1        82        82        80        80

as.data.frame(
  lapply(list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4)), function(i) {
    zz <- GenomicRanges::intersect(z[[i[1]]], z[[i[2]]])
    val <- list(round(sum(width(zz)) / sum(width(z[[i[1]]])) * 100, 0),
                round(sum(width(zz)) / sum(width(z[[i[2]]])) * 100, 0))
    setNames(val, names(z)[c(i[1], i[2])])
  }))
#>                            AH46921 AH46922 AH46921.1 AH46924 AH46921.2
#> Active TSS                      66      82        65      80        71
#> Bivalent Enhancer               28      14        32      12        22
#> Bivalent/Poised TSS             46      28        39      24        45
#> Enhancers                       57      58        62      53        42
#> Flanking Active TSS             52      59        59      52        42
#> Flanking Bivalent TSS/Enh       42      25        59      18        30
#> Genic enhancers                 41      44        46      39        35
#> Heterochromatin                 36      63        49      54        32
#> Quiescent/Low                   89      90        88      91        92
#> Repressed PolyComb              73      35        64      35        60
#> Strong transcription            65      79        71      74        58
#> Transcr. at gene 5' and 3'      39      48        45      40        38
#> Weak Repressed PolyComb         59      38        58      40        51
#> Weak transcription              67      71        67      70        66
#> ZNF genes & repeats             56      73        69      60        49
#>                            AH46926 AH46922.1 AH46924.1 AH46922.2 AH46926.1
#> Active TSS                      77        76        76        81        71
#> Bivalent Enhancer               15        31        22        20        26
#> Bivalent/Poised TSS             32        45        45        43        50
#> Enhancers                       59        68        57        44        61
#> Flanking Active TSS             62        65        51        46        60
#> Flanking Bivalent TSS/Enh       27        58        30        27        40
#> Genic enhancers                 41        54        43        38        42
#> Heterochromatin                 61        58        37        43        48
#> Quiescent/Low                   89        90        92        93        89
#> Repressed PolyComb              46        53        60        42        68
#> Strong transcription            81        79        67        66        75
#> Transcr. at gene 5' and 3'      43        51        38        42        40
#> Weak Repressed PolyComb         44        57        60        46        61
#> Weak transcription              69        71        71        70        69
#> ZNF genes & repeats             75        77        51        60        69
#>                            AH46924.2 AH46926.2
#> Active TSS                        77        68
#> Bivalent Enhancer                 15        28
#> Bivalent/Poised TSS               36        42
#> Enhancers                         37        62
#> Flanking Active TSS               39        64
#> Flanking Bivalent TSS/Enh         20        56
#> Genic enhancers                   32        44
#> Heterochromatin                   34        57
#> Quiescent/Low                     93        87
#> Repressed PolyComb                43        60
#> Strong transcription              59        79
#> Transcr. at gene 5' and 3'        33        43
#> Weak Repressed PolyComb           48        60
#> Weak transcription                66        65
#> ZNF genes & repeats               45        77
```

Consequently, it is useful to stratify our DMRs and DAPs by the directionality of the `meanDiff` and `logFC` when comparing to chromHMM tracks.

-   hypo-DMRs: Those DMRs that are hypomethylated in NA\_pos relative to BA9\_pos
-   hyper-DMRs: Those DMRs that are hypermethylated in NA\_pos relative to BA9\_pos

There are more hyper-DMRs (n = 9584) than hypo-DMRs (n = 3311). Similarly, there are more hyper-DMR-CpGs (n = 198955) than hypo-DMR-CpGs (n = 52269).

**NOTE:** An *UpSet* plot isn't perhaps appropriate/necessary since chromHMM categories are disjoint by definition and each CpG can only fall in a single category. Could just take upper bar plot for figure.

We run two analyses:

1.  Using the 13,074 POS DMRs
2.  Using the 12,895 NA\_pos vs. BA9\_pos DMRs accounting for the direction of the DMR

DMR-CpGs
--------

### POS DMRs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-2-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-2-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-2-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-2-4.png)

### Hypo NA\_pos vs. BA9\_pos DMRs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-3-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-3-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-3-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-3-4.png)

### Hyper NA\_pos vs. BA9\_pos DMRs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-4-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-4-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-4-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-4-4.png)

bigDMR-CpGs
-----------

### POS DMRs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-5-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-5-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-5-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-5-4.png)

### Hypo NA\_pos vs. BA9\_pos DMRs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-6-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-6-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-6-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-6-4.png)

### Hyper NA\_pos vs. BA9\_pos DMRs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-7-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-7-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-7-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-7-4.png)

Summary
-------

### POS DMRs

-   Choice of cell type doesn't make a big difference for top 4 categories
-   Most DMR-CpGs and bigDMR-CpGs hit 'Quiescent/Low' followed by 'Weak Transcription', 'Enhancers' and 'Weak Repressed PolyComb'
-   Not totally surprising given that 'Quiescent/Low' and 'Weak transcription' are the two largest states
-   'Enhancers' are more interesting

### NA\_pos vs. BA9\_pos DMRs

-   Choice of cell type doesn't make a big difference for top 4-5 categories
-   Most hypo-DMR-CpGs, hyper-DMR-CpGs, hypo-bigDMR-CpGs, and hyper-bigDMR-CpGs hit 'Quiescent/Low' followed by 'Weak Transcription' and 'Enhancers'
-   Not totally surprising given that 'Quiescent/Low' and 'Weak transcription' are the two largest states
-   'Enhancers' are more interesting

chromHMM context of ATAC peaks, DAPs, and bigDAPs
=================================================

Like for NA\_pos vs. BA9\_pos DMRs, we stratify DAPs and bigDAPs by hypo and hyper. It is not necessary to stratify the overall set of ATAC-seq peaks.

ATAC-peaks
----------

![](chromHMM_files/figure-markdown_github/unnamed-chunk-8-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-8-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-8-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-8-4.png)

DAPs
----

### Hypo NA\_pos vs. BA9\_pos DAPs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-9-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-9-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-9-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-9-4.png)

### Hyper NA\_pos vs. BA9\_pos DAPs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-10-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-10-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-10-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-10-4.png)

bigDAPs
-------

### Hypo NA\_pos vs. BA9\_pos DAPs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-11-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-11-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-11-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-11-4.png)

### Hyper NA\_pos vs. BA9\_pos DAPs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-12-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-12-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-12-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-12-4.png)

Summary
-------

-   Choice of cell type doesn't make a big difference (perhaps a little bigger than other in above)
-   Most hypo-DAPs, hyper-DAPs, hypo-bigDAPs, and hyper-bigDAPs hit 'Quiescent/Low' followed by 'Weak Transcription', 'Enhancers', and 'Weak Repressed PolyComb'
-   Not totally surprising given that 'Quiescent/Low', 'Weak transcription', and 'Weak Repressed Polycomb' are the two largest states
-   'Enhancers' are a little more interesting

chromHMM enrichment/depletion of DMR-CpGs and bigDMR-CpGs
=========================================================

POS DMRs
--------

![](chromHMM_files/figure-markdown_github/unnamed-chunk-13-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-13-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-13-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-13-4.png)

NA\_pos vs. BA9\_pos DMRs
-------------------------

![](chromHMM_files/figure-markdown_github/unnamed-chunk-14-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-14-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-14-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-14-4.png)

Summary
-------

### POS DMRs

-   DMR-CpGs vs. bigDMR-CpGs

### NA\_pos vs. BA9\_pos DMRs

-   DMR-CpGs vs. bigDMR-CpGs
-   Direction of signal is almost always the same
-   Magnitude of signal is also fairly similar
-   **Removed bigDMR-CpGs from plot to better see DMR-CpGs**
-   Choice of chromHMM track
-   All four tracks give very similar results, therefore reasonable to average over these
-   Based on average track
-   Most enriched: 'Transc. at a gene 5' and 3'', 'Genic enhancers', 'Flanking Active TSS', 'Enhancers', 'Bivalent Enhancers'
-   Most depleted: 'Heterochromatin', 'ZNF genes & repeats', 'Quiescent/Low'

chromHMM enrichment/depletion of ATAC-seq peaks, DAPs, and bigDAPs
==================================================================

We consider enrichment of:

1.  Peaks vs. rest of genome
2.  DAPS vs. rest of genome
    1.  bigDAPs vs. rest of genome

3.  DAPs vs. null-peaks
    1.  bigDAPs vs. null-peaks

4.  DAPs vs. non-DAPs
    1.  bigDAPs vs. non-DAPs

Note that the non-DAPs will still have some 'differential' peaks (`adj.P.Val < 0.05` but with a `abs(logFC) < 1`) whereas the null-peaks are those with `adj.P.Val > 0.05`.

We consider two different ways to calculate enrichment using ATAC-seq data:

A. Counting peaks (only 3, 4) B. Counting bases (1, 2, 3, 4)

Like for NA\_pos vs. BA9\_pos DMRs, we stratify DAPs and bigDAPs by hypo and hyper.

Counting peaks
--------------

### DAPs vs. null-peaks

![](chromHMM_files/figure-markdown_github/unnamed-chunk-15-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-15-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-15-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-15-4.png)

### DAPs vs. non-DAPs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-16-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-16-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-16-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-16-4.png)

Counting bases
--------------

### Peaks vs. rest of genome

![](chromHMM_files/figure-markdown_github/unnamed-chunk-17-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-17-2.png)

### DAPs vs. rest of genome

![](chromHMM_files/figure-markdown_github/unnamed-chunk-18-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-18-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-18-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-18-4.png)

### DAPs vs. null-peaks

![](chromHMM_files/figure-markdown_github/unnamed-chunk-19-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-19-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-19-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-19-4.png)

### DAPs vs. non-DAPs

![](chromHMM_files/figure-markdown_github/unnamed-chunk-20-1.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-20-2.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-20-3.png)![](chromHMM_files/figure-markdown_github/unnamed-chunk-20-4.png)
