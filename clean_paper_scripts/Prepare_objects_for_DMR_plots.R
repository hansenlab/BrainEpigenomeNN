# Prepare objects for plots of DMRs to be used in paper
# Only data within 1 MB of a DMR or block is retained since this is sufficient
# for plotting purposes
# Peter Hickey
# 2017-02-03

###-----------------------------------------------------------------------------
### Load packages
###

library(GenomicRanges)
library(bsseq)
library(dplyr)
library(purrr)

###-----------------------------------------------------------------------------
### Pseudocode from Lindsay for obtaining the DMRs and blocks to be used in
### plots from paper (/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/DMRs_BLOCKs_for_plotting.rda)
###

# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/Sig_F-stat_DMR_Objects.rda")
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/RefSeq_exons3.rds") # will likely change this fild
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/Sig_F-statDMRs_with_geneNames.rda")  # This also has ATAC_enh track for annoTrack plotting
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/gr.enhs.FANTOM5.rda")
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/Brain_chromHMM_Features_for_OR_calc.rda")
#
# ### load BLOCK POS DMRS ###
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/All_BLOCK_POS_DMRs_fwer50.rda")
# sig_pos_block_dmrs_gr=sig_block_dmrs_gr
#
# #### LOAD ALL BLOCK DMRs
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/All_BLOCK_DMRs_fwer50.rda")
# sig_block_dmrs_gr=GRanges(sig_block_dmrs)
#
# ### Load nonNA POS DMRs ###
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/All_NeuNpos_nonNA_DMRs_fwer50.rda")
#
# ## Load BS objects
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.small.somatic.all.pos.rda")  #Just POS
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.large.somatic.all.rda")      # for BLOCKs
# load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/BS.fit.small.sorted.somatic.all.rda")  #all samples
# #### LOAD BRAIN SPECIFIC ENHANCERS ####
#
# Brain_enhancers2=read.csv("/dcl01/feinberg/data/personal/lrizzard/ENCODE_enhancers/CelRep_enhancers.csv",header=FALSE,sep=",",skip=1)
# colnames(Brain_enhancers2)=c("chr","start","end")
# Brain_enhancers=GRanges(Brain_enhancers2) #83553
# #Brain_enhancers3=dropSeqlevels(Brain_enhancers,c("chrY","chrX","chrM"))  #used this one for OR calculation
#
#
# ### LOAD ATAC PEAKS #####
# ATAC=read.csv("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/topTable.NA_posvsBA9_pos.ATAC-seq.csv")
# colnames(ATAC)=c("X","chr","start","end","width","strand","logFC","AveExpr","t","P.Value","adj.P.Val","B")
# ATAC_gr=GRanges(ATAC)
# sig_ATAC_gr=ATAC_gr[which(ATAC_gr$adj.P.Val<=0.05)]  #71346
#
# #### LOAD DE GENES #####
#
# DE_genes=read.csv("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/topTable.NA_posvsBA9_pos.RNA-seq.csv")
# library(biomaRt)
# #mart=useDataset("hsapiens_gene_ensembl", useMart("ensembl"))  #gets GrCH38 not 37(ie hg19)
# names=strsplit(as.character(DE_genes[,1]),"[.]")
# ENSEMBL=sapply(names,"[[",1)
# DE_genes[,1]=ENSEMBL
# grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
# x=getBM(
#   filters= "ensembl_gene_id",
#   attributes= c("ensembl_gene_id", "external_gene_name", "entrezgene","description"),
#   values= as.character(DE_genes[,1]),
#   mart= grch37)
# colnames(x)=c("X","external_gene_name" ,"entrezgene","description")
# y=merge(DE_genes,x,by="X")
# sig_DE_genes=subset(y,adj.P.Val<=0.05) #3152
# sig_DE_genes2=sig_DE_genes[with(sig_DE_genes,order(adj.P.Val)),]
#
#
#
# ###### GET LIST OF DMRs to plot
# nonNA=as.data.frame(sig_pos_nonNA_dmrs_gr)
#
# LRRC4=GRanges(sig_pos_dmrs[which(sig_pos_dmrs$start==127668929),])
#
#
# #GRIK4=nonNA[which(nonNA$start==120703308),]   #extend 5000   COME BACK TO THIS ONE...not in original 208 DMR list
# # instead of GRIK4 can use this one ......
#
# DLGAP2=GRanges(nonNA[which(nonNA$start==1617411),])   # extend 5000
#
#
# GABRB2_block=as.data.frame(sig_pos_block_dmrs_gr)[which(as.data.frame(sig_pos_block_dmrs_gr)$start==160765325),]
# GABRB2_block=GRanges(GABRB2_block) #extend 100,000
# RBFOX3=GRanges(sig_dmrs[which(sig_dmrs$start==77461456),]) #extend 5000
# QKI_block=GRanges(sig_block_dmrs[which(sig_block_dmrs$start==163831857),]) #extend 100,000
# MEF2C_block=GRanges(as.data.frame(sig_pos_block_dmrs_gr)[which(as.data.frame(sig_pos_block_dmrs_gr)$start==87979404),]) # extend 100,000
# MEF2C_small_DMR=GRanges(sig_pos_dmrs[which(sig_pos_dmrs$start==88027045),]) # extend 5000 (this one is optional)
# RGS9=GRanges(sig_pos_dmrs[which(sig_pos_dmrs$start==63134131),]) # extend 5000
# SATB2=GRanges(sig_pos_dmrs[which(sig_pos_dmrs$start==200246174),]) # extend 5000
# SEMA7A=GRanges(sig_pos_dmrs[which(sig_pos_dmrs$start==74724026),]) # extend 5000
#
# DMRs_BLOCKs_for_plotting=list("LRRC4"=LRRC4,"DLGAP2"=DLGAP2,"GABRB2_block"=GABRB2_block,"RBFOX3"=RBFOX3,"QKI_block"=QKI_block,
#                               "MEF2C_block"=MEF2C_block,"MEF2C_small_DMR"=MEF2C_small_DMR,"RGS9"=RGS9,"SATB2"=SATB2,"SEMA7A"=SEMA7A)
# setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/")
# save(DMRs_BLOCKs_for_plotting,file="DMRs_BLOCKs_for_plotting.rda")


###-----------------------------------------------------------------------------
### Create objects
###

#-------------------------------------------------------------------------------
# All DMRs and BLOCKS within POS
#

# DMRs
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/POS_DMRs_fwer50_CpG5.rda")
pos_dmrs <- GRanges(sig_pos_dmrs)

# Blocks
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/All_BLOCK_POS_DMRs_fwer50.rda")
pos_blocks <- GRanges(sig_block_dmrs)

#-------------------------------------------------------------------------------
# Load all POSvsNEG DMRs
#

load("../Objects/All_Annotated_DMRs_GRanges.rda")
pos_vs_neg_dmrs <- Annotated_POSvNEG_DMRs_gr

#-------------------------------------------------------------------------------
# Load all POSvsNEG BLOCKs
#

load("../Objects/Annotated_POSvNEG_BLOCKs_GRanges.rda")
pos_vs_neg_blocks <- Annotated_POSvNEG_BLOCKs_gr

#-------------------------------------------------------------------------------
# Load all non-NA, NeuN+ DMRs
#

load("../Objects/All_NeuNpos_nonNA_DMRs_fwer50.rda")
pos_non_NA_dmrs <- sig_pos_nonNA_dmrs_gr

#-------------------------------------------------------------------------------
# Load assays-and-features.rda to get the rna_atac_meth object
#

load("../integrating-dmrs-daps-and-degs/objects/assays-and-features.rda")
rna_atac_meth <- rna_atac_meth %>%
  mutate(npD = pmap_int(list(pDAP, pDMR), sum),
         neD = pmap_int(list(eDAP, eDMR), sum),
         npeD = npD + neD)

sameSign <- function(x, i) {
  s <- sign(x[i])
  if (length(s)) {
    all(s == s[1])
  } else {
    TRUE
  }
}

# Pick some genes from rna_atac_meth with neD > 1 and are consistent
enhancer_rich_consistent_DE_genes <- rna_atac_meth %>%
  filter(DE, gene %in% names(fantom5_enhancers_by_gene), neD > 0) %>%
  mutate(ss = map2_lgl(
    pmap(list(pLogFC, eLogFC, map(map2(pMeanDiff, eMeanDiff, c), `-`),
              as.list(expLogFC)), c),
    pmap(list(pDAP, eDAP, pDMR, eDMR, as.list(DE)), c),
    sameSign)) %>%
  filter(ss) %>%
  arrange(desc(neD)) %>%
  filter(neD > 1) %>%
  .$gene

# Create a GRanges containing the gene body and all its linked enhancers
enhancer_rich_consistent_DE_genes_gr <-
  lapply(enhancer_rich_consistent_DE_genes, function(x) {
    y <- c(granges(unflattened_features$genes[x]),
           granges(fantom5_enhancers_by_gene[[x]]))
    y <- GenomicRanges::reduce(unstrand(y))
    GRanges(seqnames(y)[1],
            IRanges(min(start(y)),
                    max(end(y))))
  })

enhancer_rich_consistent_DE_genes_gr <-
  unlist(GRangesList(enhancer_rich_consistent_DE_genes_gr))
names(enhancer_rich_consistent_DE_genes_gr) <-
  enhancer_rich_consistent_DE_genes
enhancer_rich_consistent_DE_genes_gr <-
  resize(enhancer_rich_consistent_DE_genes_gr,
         width = width(enhancer_rich_consistent_DE_genes_gr) + 2 * 5000,
         fix = "center")
list_of_enhancer_rich_consistent_DE_genes_gr <-
  as.list(split(enhancer_rich_consistent_DE_genes_gr,
                names(enhancer_rich_consistent_DE_genes_gr)))

#-------------------------------------------------------------------------------
# Regions to plot (split into DMRs and blocks)
#

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/DMRs_BLOCKs_for_plotting.rda")
DMRs_blocks_for_plotting <- c(DMRs_BLOCKs_for_plotting,
                              list("SIX3" =  GRanges("chr2:45170354-45171401")),
                              list_of_enhancer_rich_consistent_DE_genes_gr)
DMRs_blocks_for_plotting <- GRangesList(endoapply(DMRs_blocks_for_plotting,
                                                  granges))
blocks <- DMRs_blocks_for_plotting[grep("block",
                                        names(DMRs_blocks_for_plotting))]
dmrs <- DMRs_blocks_for_plotting[grep("block",
                                      names(DMRs_blocks_for_plotting),
                                      invert = TRUE)]
dmrs_resized <- resize(dmrs, width(dmrs) +  2 * 1000000, fix = "center")
blocks_resized <- resize(blocks, width(blocks) +  2 * 1000000, fix = "center")

#-------------------------------------------------------------------------------
# BSseq objects
#

# Small smooth
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.small.somatic.all.pos.rda")
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.small.somatic.all.neg.rda")
BSseq_small_smooth <- cbind(BSseq.pos, BSseq.neg)
rm(BSseq.pos, BSseq.neg)

# Large smooth
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.large.somatic.all.rda")
BSseq_large_smooth <- BSseq.all
rm(BSseq.all)

BSseq_dmrs <- subsetByOverlaps(BSseq_small_smooth, dmrs_resized)
BSseq_blocks <- subsetByOverlaps(BSseq_large_smooth, blocks_resized)

#-------------------------------------------------------------------------------
# ATAC-seq coverage SummarizedExperiment objects and ATAC-seq peaks as GRanges
# objects

atac_cov <- readRDS("/dcl01/hansen/data/flow-sorted-brain-atac/objects/flow-sorted-brain-atac.cov.list_of_se.rds")
# Drop 'rep2' samples (technical replicates)
atac_cov <- lapply(atac_cov, function(se) {
  se[, se$REPLICATE == "rep1"]
})

atac_cov_dmrs <- lapply(atac_cov, function(se) {
  subsetByOverlaps(se, dmrs_resized)
})
atac_cov_blocks <- lapply(atac_cov, function(se) {
  subsetByOverlaps(se, blocks_resized)
})

load("../Objects/atac.rda")
atac_overall_peaks_dmrs <- subsetByOverlaps(atac$overall_peaks, dmrs_resized)
atac_overall_peaks_blocks <- subsetByOverlaps(atac$overall_peaks,
                                              blocks_resized)
atac_NA_posvsBA9_pos_dmrs <- subsetByOverlaps(atac$NA_posvsBA9_pos,
                                              dmrs_resized)
atac_NA_posvsBA9_pos_blocks <- subsetByOverlaps(atac$NA_posvsBA9_pos,
                                                blocks_resized)
atac_NA_negvsBA9_neg_dmrs <- subsetByOverlaps(atac$NA_negvsBA9_neg,
                                              dmrs_resized)
atac_NA_negvsBA9_neg_blocks <- subsetByOverlaps(atac$NA_negvsBA9_neg,
                                                blocks_resized)
atac_ave_pos_vs_ave_neg_dmrs <- subsetByOverlaps(atac$ave_pos_vs_ave_neg,
                                                 dmrs_resized)
atac_ave_pos_vs_ave_neg_blocks <- subsetByOverlaps(atac$ave_pos_vs_ave_neg,
                                                   blocks_resized)

# Only retain those DAPs with adjusted P-value < 0.05
# NOTE: None of the NA_negvsBA9_neg DA peaks that are near DMRs and BLOCKs
#       have an adjusted P-value < 0.05
atac_NA_posvsBA9_pos_dmrs <-
  atac_NA_posvsBA9_pos_dmrs[atac_NA_posvsBA9_pos_dmrs$adj.P.Val < 0.05]
atac_NA_posvsBA9_pos_blocks <-
  atac_NA_posvsBA9_pos_blocks[atac_NA_posvsBA9_pos_blocks$adj.P.Val < 0.05]
atac_NA_negvsBA9_neg_dmrs <-
  atac_NA_negvsBA9_neg_dmrs[atac_NA_negvsBA9_neg_dmrs$adj.P.Val < 0.05]
atac_NA_negvsBA9_neg_blocks <-
  atac_NA_negvsBA9_neg_blocks[atac_NA_negvsBA9_neg_blocks$adj.P.Val < 0.05]
atac_ave_pos_vs_ave_neg_dmrs <-
  atac_ave_pos_vs_ave_neg_dmrs[atac_ave_pos_vs_ave_neg_dmrs$adj.P.Val < 0.05]
atac_ave_pos_vs_ave_neg_blocks <-
  atac_ave_pos_vs_ave_neg_blocks[atac_ave_pos_vs_ave_neg_blocks$adj.P.Val < 0.05]

#-------------------------------------------------------------------------------
# RefSeq gene models
#

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/RefSeq_exons3.rds")
RefSeq_exons3_dmrs <- as.data.frame(
  subsetByOverlaps(
    makeGRangesFromDataFrame(RefSeq_exons3, keep.extra.columns = TRUE),
    dmrs_resized)
)
RefSeq_exons3_blocks <- as.data.frame(
  subsetByOverlaps(
    makeGRangesFromDataFrame(RefSeq_exons3, keep.extra.columns = TRUE),
    blocks_resized)
)

#-------------------------------------------------------------------------------
# Brain enhancers
#

load("../Objects/Brain_enh.rda")
Brain_enh_dmrs <- subsetByOverlaps(Brain_enh, dmrs_resized)
Brain_enh_blocks <- subsetByOverlaps(Brain_enh, blocks_resized)

###-----------------------------------------------------------------------------
### Save objects
###

save(pos_dmrs, pos_blocks,
     pos_vs_neg_dmrs, pos_vs_neg_blocks, pos_non_NA_dmrs,
     dmrs, blocks,
     BSseq_dmrs, BSseq_blocks,
     atac_cov_dmrs, atac_cov_blocks,
     atac_overall_peaks_dmrs, atac_overall_peaks_blocks,
     atac_NA_posvsBA9_pos_dmrs, atac_NA_posvsBA9_pos_blocks,
     atac_NA_negvsBA9_neg_dmrs, atac_NA_negvsBA9_neg_blocks,
     atac_ave_pos_vs_ave_neg_dmrs, atac_ave_pos_vs_ave_neg_blocks,
     RefSeq_exons3_dmrs, RefSeq_exons3_blocks,
     Brain_enh_dmrs, Brain_enh_blocks,
     file = "../Objects/objects_for_plotting_DMR_and_block_examples.rda")
