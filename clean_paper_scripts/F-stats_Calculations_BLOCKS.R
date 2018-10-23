### Perform F-stat calculations with 1000 permutations for all comparisons
# qrsh -l mem_free=70G,h_vmem=71G,h_fsize=200G -pe local 4
library(devtools)
devtools::install_github("PeteHaitch/bsseq")
library(bsseq)
devtools::session_info()

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.large.rda")
BSseq <- BS.fit.large
BSseq <- sort(sortSeqlevels(BSseq))
BSseq=dropSeqlevels(BSseq,c("chrY","chrX","chrM","lambda"))
BS.cov=getCoverage(BSseq)
keepLoci.all <- which(rowSums(BS.cov >= 1) >= 45)

BSseq.all<-BSseq[keepLoci.all,]  #23059530 loci

pryr::object_size(BSseq.all)  #25.1G
rm(BS.fit.large,BSseq)
gc()
as.data.frame(pData(BSseq.all))

pData(BSseq.all)$Tissue_NeuN <- paste(pData(BSseq.all)$Tissue, pData(BSseq.all)$NeuN,sep = "_")
design <- model.matrix(~BSseq.all$Tissue * BSseq.all$NeuN)
colnames(design) <- gsub("BSseq.all\\$", "", colnames(design))
contrasts <- diag(rep(1,8))
contrasts <- contrasts[, -1]
rownames(contrasts) <- colnames(design)
fac <- pData(BSseq.all)$Tissue_NeuN

set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.all, design = design,
                               contrasts = contrasts,
                               cutoff = 2^2, fac = fac, type = "blocks",
                               nperm = 1000, coef = NULL, maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 10000, mc.cores = 6)

setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/")

save(dmrs,file="F-stat_1000perm_BLOCK_dmrs.rda")

dmr_list=dmrs$dmrs # DMRs
bstat=dmrs$bstat
sig_block_dmrs=dmr_list[which(dmr_list$fwer<=50),] #20373
sig_block_dmrs_gr=GRanges(sig_block_dmrs)

save(sig_block_dmrs,file="All_BLOCK_DMRs_fwer50.rda")

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/RefSeq_exons3.rds")
gt=RefSeq_exons3
pData<-pData(BSseq.all)

pData$Sex_color=pData$Sex
pData$Sex_color=replace(pData$Sex_color,pData$Sex_color=="female","firebrick3")
pData$Sex_color=replace(pData$Sex_color,pData$Sex_color=="male","deepskyblue")
pData$Tissue_color=pData$Tissue
pData$Tissue_color=replace(pData$Tissue_color,pData$Tissue_color=="NA","chocolate1")
pData$Tissue_color=replace(pData$Tissue_color,pData$Tissue_color=="BA9","deepskyblue")
pData$Tissue_color=replace(pData$Tissue_color,pData$Tissue_color=="BA24","deeppink")
pData$Tissue_color=replace(pData$Tissue_color,pData$Tissue_color=="HC","darkgrey")
pData$lty<- ifelse(pData$NeuN == "pos", 1, 2)
pData(BSseq.all)=pData


pdf("BLOCK_DMRs_nperm1000_fwer50.pdf",width=10,height=5)
bsseq:::plotManyRegions(BSseq.all,sig_block_dmrs_gr[1:100,],extend=20000,addRegions=sig_block_dmrs_gr,col=pData(BSseq.all)$Tissue_color,geneTrack=gt,main="BLOCK_DMRs")
dev.off()

sig_block_dmrs_10percent=sig_block_dmrs_gr[which(sig_block_dmrs_gr$maxDiff>=0.1),] #4185
pdf("BLOCK_DMRs_10percentDiff_fwer50.pdf",width=10,height=5)
bsseq:::plotManyRegions(BSseq.all,sig_block_dmrs_10percent[1:100,],extend=20000,addRegions=sig_block_dmrs_10percent,col=pData(BSseq.all)$Tissue_color,geneTrack=gt,main="BLOCK_DMRs >=10% Diff")
dev.off()


########### DMR BLOCKS with just POS ###############
### DO WHEN CLUSTER BACK ONLINE ######

names=sampleNames(BSseq.all)
pos=names[grep("pos",names)]
BSseq.pos=BSseq.all[,pos]


design <- model.matrix(~BSseq.pos$Tissue)
colnames(design) <- gsub("BSseq.pos\\$Tissue", "", colnames(design))

contrasts <- cbind(c(0,1,0,0),
                   c(0,0,1,0),
                   c(0,0,0,1))

rownames(contrasts) <- colnames(design)
fac <- pData(BSseq.pos)$Tissue

set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.pos, design = design,
                               contrasts = contrasts,
                               cutoff = 2 ^ 2, fac = fac,
                               nperm = 1000, coef = NULL, maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 10000, mc.cores = 6)


setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/")

save(dmrs,file="F-stat_1000perm_BLOCK_POS_dmrs.rda")

dmr_list=dmrs$dmrs # DMRs
bstat=dmrs$bstat
sig_block_dmrs=dmr_list[which(dmr_list$fwer<=50),] #1808
sig_block_dmrs_gr=GRanges(sig_block_dmrs)

save(sig_block_dmrs,file="All_BLOCK_POS_DMRs_fwer50.rda")
