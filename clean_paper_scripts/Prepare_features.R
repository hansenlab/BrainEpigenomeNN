## This script takes various genomic features
## and prepares them for overlap.
## This includes only looking at autosomes etc.

## Set of DEs (gene-level and transcript-level)
## Set of D-ATAC

## GENCODE v19: genic, exonic, intronic, promoters, 5' UTR, 3' UTR, intergenic
## CGI: islands, shores, shelves, open sea
## Brain Enhancers
## Fantom5 enhancers (TSS-linked enhancers)
## ChromHMM: 15 state model

## Published DMRs:
## 1. Carolina's paper (CHARM neuron vs glia)
## 2. Kozlenkov et al. 2014 NAR (450k neuron vs glia)
## 3. Lister Science 2013 (WGBS neuron vs glia)

library(GenomicRanges)
library(bsseq)
library(readr)
library(dplyr)
library(tibble)
autosomes <- paste0("chr", 1:22)

## GENCODE v19
## Prepared with https://github.com/feinberglabepigenetics/GTExScripts/blob/master/FlowSortingProject/genomic-features/scripts/gene-models.R
## These objects have already been restricted to the autosomes
load("../genomic-features/objects/flattened-GENCODE-v19-features.rda")
genes <- readRDS("../genomic-features/objects/GENCODE-v19-genes.rds")
transcripts_by_gene <- readRDS("../genomic-features/objects/GENCODE-v19-transcripts-by-gene.rds")
gencode_features <- list("union" = flattened_features,
                         "pc_transcripts" = flattened_features_pc_transcripts,
                         "lnc_transcripts" = flattened_features_lnc_transcripts,
                         "genes" = genes,
                         "transcripts_by_gene" = transcripts_by_gene)
save(gencode_features,
     file = "../Objects/gencode_features.rda")

## Set of DEs

# NA_pos vs BA9_pos
gene_level_de_pos <-
  read_csv("../RNA-seq/extdata/topTable.NA_posvsBA9_pos.RNA-seq.csv.gz")
colnames(gene_level_de_pos)[1] <- "gene_id"
gene_level_de_pos_gr <- gene_level_de_pos %>%
  inner_join(rownames_to_column(as.data.frame(gencode_features$genes)),
             by = c("gene_id" = "rowname")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
transcript_level_de_pos <-
  read_csv("../RNA-seq/extdata/topTable.NA_posvsBA9_pos.RNA-seq-transcripts.csv.gz")
colnames(transcript_level_de_pos)[1] <- "tx_name"
transcript_level_de_pos_gr <- transcript_level_de_pos %>%
  inner_join(as.data.frame(unlist(
    gencode_features$transcripts_by_gene, use.names = FALSE)),
             by = c("tx_name" = "tx_name")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
rna_seq_de_pos <- list("gene_level" = gene_level_de_pos_gr,
               "transcript_level" = transcript_level_de_pos_gr)

# NA_neg vs BA9_neg
gene_level_de_neg <-
  read_csv("../RNA-seq/extdata/topTable.NA_negvsBA9_neg.RNA-seq.csv.gz")
colnames(gene_level_de_neg)[1] <- "gene_id"
gene_level_de_neg_gr <- gene_level_de_neg %>%
  inner_join(rownames_to_column(as.data.frame(gencode_features$genes)),
             by = c("gene_id" = "rowname")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
transcript_level_de_neg <-
  read_csv("../RNA-seq/extdata/topTable.NA_negvsBA9_neg.RNA-seq-transcripts.csv.gz")
colnames(transcript_level_de_neg)[1] <- "tx_name"
transcript_level_de_neg_gr <- transcript_level_de_neg %>%
  inner_join(as.data.frame(unlist(
    gencode_features$transcripts_by_gene, use.names = FALSE)),
    by = c("tx_name" = "tx_name")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
rna_seq_de_neg <- list("gene_level" = gene_level_de_neg_gr,
                       "transcript_level" = transcript_level_de_neg_gr)

# ave_pos vs ave_neg
gene_level_de_ave_pos_vs_ave_neg <-
  read_csv("../RNA-seq/extdata/topTable.ave_pos_vs_ave_neg.RNA-seq.csv.gz")
colnames(gene_level_de_ave_pos_vs_ave_neg)[1] <- "gene_id"
gene_level_de_ave_pos_vs_ave_neg_gr <- gene_level_de_ave_pos_vs_ave_neg %>%
  inner_join(rownames_to_column(as.data.frame(gencode_features$genes)),
             by = c("gene_id" = "rowname")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
transcript_level_de_ave_pos_vs_ave_neg <-
  read_csv("../RNA-seq/extdata/topTable.ave_pos_vs_ave_neg.RNA-seq-transcripts.csv.gz")
colnames(transcript_level_de_ave_pos_vs_ave_neg)[1] <- "tx_name"
transcript_level_de_ave_pos_vs_ave_neg_gr <- transcript_level_de_ave_pos_vs_ave_neg %>%
  inner_join(as.data.frame(unlist(
    gencode_features$transcripts_by_gene, use.names = FALSE)),
    by = c("tx_name" = "tx_name")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
rna_seq_de_ave_pos_vs_ave_neg <-
  list("gene_level" = gene_level_de_ave_pos_vs_ave_neg_gr,
       "transcript_level" = transcript_level_de_ave_pos_vs_ave_neg_gr)

save(rna_seq_de_pos,
     rna_seq_de_neg,
     rna_seq_de_ave_pos_vs_ave_neg,
     file = "../Objects/rna_seq.rda")

## Set of ATAC-seq peaks and differentially accessibly peaks (NA_posvsBA9_pos,
## NA_negvsBA9_neg, ave_pos_vs_ave_neg)
atac_overall_peaks_gr <- keepSeqlevels(
  rowRanges(readRDS("../ATAC-seq/extdata/flow-sorted-brain-atac/objects/flow-sorted-brain-atac.overall.se.rds")),
  autosomes)
atac_NA_posvsBA9_pos <-
  read_csv("../ATAC-seq/extdata/topTable.NA_posvsBA9_pos.ATAC-seq.csv.gz")
atac_NA_negvsBA9_neg <-
  read_csv("../ATAC-seq/extdata/topTable.NA_negvsBA9_neg.ATAC-seq.csv.gz")
atac_ave_pos_vs_ave_neg <-
  read_csv("../ATAC-seq/extdata/topTable.ave_pos_vs_ave_neg.ATAC-seq.csv.gz")

atac_NA_posvsBA9_pos_gr <- GRanges(atac_NA_posvsBA9_pos[, -1])
atac_NA_negvsBA9_neg_gr <- GRanges(atac_NA_negvsBA9_neg[, -1])
atac_ave_pos_vs_ave_neg_gr <- GRanges(atac_ave_pos_vs_ave_neg[, -1])

atac <- list("overall_peaks" = atac_overall_peaks_gr,
             "NA_posvsBA9_pos" = atac_NA_posvsBA9_pos_gr,
             "NA_negvsBA9_neg" = atac_NA_negvsBA9_neg_gr,
             "ave_pos_vs_ave_neg" = atac_ave_pos_vs_ave_neg_gr)
save(atac,
     file = "../Objects/atac.rda")

## CGI
## From UCSC and prepared with
## https://github.com/feinberglabepigenetics/GTExScripts/blob/master/FlowSortingProject/genomic-features/scripts/cgi-related-features.R)
load("../genomic-features/objects/CGI-related-features-hg19.rda")
CGI <- keepSeqlevels(cgis, autosomes)
shelves <- keepSeqlevels(shelves, autosomes)
shores <- keepSeqlevels(shores, autosomes)
open_sea <- keepSeqlevels(open_sea, autosomes)

cgi_features <- list("CGI" = CGI,
                     "Shores" = shores,
                     "Shelves" = shelves,
                     "OpenSea" = open_sea)

save(cgi_features,
     file = "../Objects/cgi_features.rda")

## Brain Enhancers
## From Vermunt 2014 Cell Reports
## H3K27ac in 136 brain regions, unsorted adult human brain

brainDF <- read.csv("/dcl01/feinberg/data/personal/lrizzard/ENCODE_enhancers/CelRep_enhancers.csv",
                    header = FALSE, sep = ",", skip = 1,
                    stringsAsFactors = FALSE)
colnames(brainDF) <- c("chr","start","end")
Brain_enh <- data.frame2GRanges(brainDF, keepColumns = TRUE) #83553
Brain_enh <- keepSeqlevels(Brain_enh, autosomes)
save(Brain_enh, file = "../Objects/Brain_enh.rda")

## (FANTOM5 TSS-enhancer associated enhancer, GENCODE v19 gene ID)-pairs)
FANTOM5_enh <-
  readRDS("../genomic-features/objects/FANTOM5_enhancers_with_GENCODE_v19_links.rds")

FANTOM5_enh <- keepSeqlevels(FANTOM5_enh, autosomes)
save(FANTOM5_enh, file = "../Objects/FANTOM5_enh.rda")

## ChromHMM

library(AnnotationHub)
ah <- AnnotationHub()

# Query AnnotationHub for all brain chromHMM tracks
queries <- query(ah, c("chromHMM", "brain"))   # 15-state model
queries <- queries[c(1:3,5:8)]
## E071 Brain Hippocampus Middle
## E074 Brain Substantia Nigra
## E068 Brain Anterior Caudate
## E069 Brain Cingulate Gyrus
## E072 Brain Inferior Temporal Lobe
## E067 Brain Angular Gyrus
## E073 Brain_Dorsolateral_Prefrontal_Cortex
## AH46920 | E067_15_coreMarks_mnemonics.bed.gz
## AH46921 | E068_15_coreMarks_mnemonics.bed.gz
## AH46922 | E069_15_coreMarks_mnemonics.bed.gz
## AH46924 | E071_15_coreMarks_mnemonics.bed.gz
## AH46925 | E072_15_coreMarks_mnemonics.bed.gz
## AH46926 | E073_15_coreMarks_mnemonics.bed.gz
chromHMM_info <- matrix(c(
  "AH46920", "E067", "AngularGyrus", "",
  "AH46921", "E068", "AnteriorCaudate",  "Adjacent to Nacc",
  "AH46922", "E069", "CingulateGyrus", "BA24 is a subset of this region",
  "AH46924", "E071", "HippocampusMiddle", "HC",
  "AH46925", "E072", "InferiorTemporalLobe", "",
  "AH46926", "E073", "DorsolateralPrefrontalCortex", "BA9",
  "AH46927", "E074", "SubstantiaNigra", ""), ncol = 4, byrow = TRUE)
chromHMM_info <- as.data.frame(chromHMM_info, stringsAsFactors = FALSE)
colnames(chromHMM_info) <- c("AH", "E", "Region", "Info")

## Download data
chromHMM_tracks <- lapply(queries, function(x) {
  x[[1]]
})
chromHMM_tracks <- lapply(chromHMM_tracks, function(x){
    keepSeqlevels(x, autosomes)
})
save(chromHMM_tracks, chromHMM_info, file = "../Objects/chromHMM.rda")


## ## For each chromHMM track, split by annotation type
## split_chromHMM_tracks <- lapply(chromHMM_tracks, function(x) {as.list(split(x, x$abbr))})
## ## Make into a single list with name of the form
## ## <AnnotationHubID>.<chromHMMState>
## split_chromHMM_tracks <- unlist(split_chromHMM_tracks, recursive = FALSE)
## split_chromHMM_tracks_autosomes<-lapply(split_chromHMM_tracks,function(x){keepSeqlevels(x, autosomes))})

## ###### Now calculate OR over ChromHMM states for each tissue type and each DMR list
## chromHMM_features = split_chromHMM_tracks_autosomes
## sub("", ""names(chromHMM_features))

