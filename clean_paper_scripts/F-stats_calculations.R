### Perform F-stat calculations with 1000 permutations for all comparisons


library(minfi)
library(devtools)
install_github("kasperdanielhansen/bsseq")
library(bsseq)
devtools::session_info()

file.copy(from = "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/BS.fit.small.sorted.somatic.all.rda",
          to = "/scratch/temp/BS.fit.small.sorted.somatic.all.rda")
load("/scratch/temp/BS.fit.small.sorted.somatic.all.rda")

unlink("/scratch/temp/BS.fit.small.sorted.somatic.all.rda")

BSseq <- BS.fit.small.sorted.somatic.all
pryr::object_size(BSseq)  #25.1G
rm(BS.fit.small.sorted.somatic.all)
gc()
as.data.frame(pData(BSseq))


#########  ALL COMPARISONS  ############

pData(BSseq)$Tissue_NeuN <- paste(pData(BSseq)$Tissue, pData(BSseq)$NeuN,sep = "_")
design <- model.matrix(~BSseq$Tissue * BSseq$NeuN)
colnames(design) <- gsub("BSseq\\$", "", colnames(design))
contrasts <- diag(rep(1,8))
contrasts <- contrasts[, -1]
rownames(contrasts) <- colnames(design)
fac <- pData(BSseq)$Tissue_NeuN

set.seed(12345)
dmrs <- fstat.pipeline(BSseq = BSseq, design = design,
                       contrasts = contrasts,
                       cutoff = c(-4.6 ^ 2, 4.6 ^ 2), fac = fac,
                       nperm = 1000, coef = NULL, maxGap.sd = 10 ^ 8,
                       maxGap.dmr = 300, mc.cores = 4)

setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/")

save(dmrs,file="f-stat_1000perm_dmrs.rda")

dmr_list=dmrs$dmrs #376711 DMRs
bstat=dmrs$bstat
sig_dmrs=dmr_list[which(dmr_list$fwer<=50),] #100875
sig_dmrs_gr=GRanges(sig_dmrs)

save(sig_dmrs,file="All_DMRs_fwer50.rda")

### DMRs with min of 5 CpGs and FWER <=50; did not use this going forward
sig_dmrs2=dmr_list[which(dmr_list$fwer<=50 & dmr_list$n>=5),] #96369

### Look at plots of best and worst DMRs based on fwer ###

pdf(file = "fwer.fstat-zero.pdf")
plotManyRegions(BSseq = BSseq,
                BSseqStat = bstat,
                regions = head(subset(dmr_list, fwer == 0), n = 100),col=pData(BSseq)$Tissue_color,
                stat = "stat", addRegions = dmr_list, extend = 10000,
                stat.ylim = c(-10, 10))
dev.off()

pdf(file = "fwer.fstat-large.pdf")
plotManyRegions(BSseq = BSseq,
                BSseqStat = bstat,
                regions = head(dmr_list[order(dmr_list$fwer, decreasing = TRUE), ],
                               n = 100),col=pData(BSseq)$Tissue_color,
                stat = "stat", addRegions = dmr_list, extend = 10000,
                stat.ylim = c(-10, 10))
dev.off()


######### JUST NEUN POSITIVE COMPARISONS  ############

names=sampleNames(BSseq)
pos=names[grepl("pos",names)] #subset just the pos
neg=names[grepl("neg",names)]

BSseq.pos=BSseq[,pos]
BSseq.neg=BSseq[,neg]
pData=pData(BSseq)


### POS first ###

design <- model.matrix(~BSseq.pos$Tissue)
colnames(design) <- gsub("BSseq.pos\\$", "", colnames(design))

contrasts <- cbind(c(0,1,0,0),
                   c(0,0,1,0),
                   c(0,0,0,1))

rownames(contrasts) <- colnames(design)
fac <- pData(BSseq.pos)$Tissue

set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.pos, design = design,
                               contrasts = contrasts,
                               cutoff = c(-4.6 ^ 2, 4.6 ^ 2), fac = fac,
                               nperm = 1000, coef = NULL, maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300, mc.cores = 4)

save(dmrs,file="POS.f-stat_1000perm_dmrs.rda")
pos_dmrs=dmrs$dmrs
sig_pos_dmrs=pos_dmrs[which(pos_dmrs$fwer<=50),]
sig_pos_dmrs_gr=GRanges(sig_pos_dmrs)

######### JUST NEUN NEGATIVE COMPARISONS  ############


design <- model.matrix(~BSseq.neg$Tissue)
colnames(design) <- gsub("BSseq.neg\\$", "", colnames(design))

contrasts <- cbind(c(0,1,0,0),
                   c(0,0,1,0),
                   c(0,0,0,1))

rownames(contrasts) <- colnames(design)

fac <- pData(BSseq.neg)$Tissue

set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.neg, design = design,
                               contrasts = contrasts,
                               cutoff = c(-4.6 ^ 2, 4.6 ^ 2), fac = fac,
                               nperm = 1000, coef = NULL, maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300, mc.cores = 4)

setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/")
neg_dmrs=dmrs
save(neg_dmrs,file="NEG.f-stat_1000perm_dmrs.rda")
nrow(neg_dmrs) #2808

sig_neg_dmrs=neg_dmrs$dmrs[which(neg_dmrs$dmrs$fwer<=50),]
sig_neg_dmrs_gr=GRanges(sig_neg_dmrs)

nrow(subset(neg_dmrs$dmrs, fwer <= 50)) # 114 DMRs
nrow(subset(neg_dmrs$dmrs, fwer <= 100)) # 147 DMRs

########## JUST UNSORTED COMPARISONS ##########
file.copy(from = "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/notsorted/BS.unsorted.fit.small.somatic.all.rda",
          to = "/scratch/temp/BS.unsorted.fit.small.somatic.all.rda")
load("/scratch/temp/BS.unsorted.fit.small.somatic.all.rda")

unlink("/scratch/temp/BS.unsorted.fit.small.somatic.all.rda")

BSseq <- BS.unsorted.fit.small.somatic.all
rm(BS.unsorted.fit.small.somatic.all)
as.data.frame(pData(BSseq))

###Leave out caudate from the analysis at this time....

BSseq2=BSseq[,which(pData(BSseq)$Tissue!="caudate")]

design <- model.matrix(~BSseq2$Tissue)
colnames(design) <- gsub("BSseq2\\$", "", colnames(design))
contrasts <- cbind(c(0,1,0,0),
                   c(0,0,1,0),
                   c(0,0,0,1))
rownames(contrasts) <- colnames(design)

fac <- pData(BSseq2)$Tissue
set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq2, design = design,
                               contrasts = contrasts,
                               cutoff = c(-4.6 ^ 2, 4.6 ^ 2), fac = fac,
                               nperm = 1000, coef = NULL, maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300, mc.cores = 4)

setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/notsorted/")
unsorted_dmrs=dmrs

save(unsorted_dmrs, file = "Unsorted_F-stat_DMRs_nperm1000.rda")

dmr_list=unsorted_dmrs$dmrs  #3059 total dmrs
sig_unsorted_dmrs=dmr_list[dmr_list$fwer<=50,] #71
sig_unsorted_dmrs_gr=GRanges(sig_unsorted_dmrs)

sig_dmrs2=dmr_list[dmr_list$fwer<=100,]  #94


###### SAVE ALL SIG DMR LISTS IN ONE OBJECT  ###########

save(sig_unsorted_dmrs_gr,sig_unsorted_dmrs,sig_pos_dmrs_gr,sig_pos_dmrs,sig_neg_dmrs_gr,sig_neg_dmrs,sig_dmrs,sig_dmrs_gr,file="Sig_F-stat_DMR_Objects.rda",compress=TRUE)

###### RESULTS FROM THIS ######

# sig_unsorted_dmrs = 71   minimum CpGs = 8
# sig_pos_dmrs = 13074     minimum CpGs = 5
# sig_neg_dmrs = 114       minimum CpGs = 11
# sig_dmrs = 100875        minimum CpGs = 2



##### SEX ANALYSIS ######
setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects")
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/BS.fit.small.sorted.somatic.all.rda")
BSseq <- BS.fit.small.sorted.somatic.all
rm(BS.fit.small.sorted.somatic.all)
### Subset by Tissue and run F-stat permutations on each tissue separately ###

BSseq_NA=BSseq[,which(pData(BSseq)$Tissue=="NA")]
BSseq_BA9=BSseq[,which(pData(BSseq)$Tissue=="BA9")]
BSseq_BA24=BSseq[,which(pData(BSseq)$Tissue=="BA24")]
BSseq_HC=BSseq[,which(pData(BSseq)$Tissue=="HC")]

save(BSseq_NA,file="BSseq_NA.rda",compress=TRUE)
save(BSseq_BA9,file="BSseq_BA9.rda",compress=TRUE)
save(BSseq_BA24,file="BSseq_BA24.rda",compress=TRUE)
save(BSseq_HC,file="BSseq_HC.rda",compress=TRUE)


#### Write this script as separate file and run for each tissue Sex_Fstat.R####
#had to screen even though only doing one tissue at a time
#qrsh -l mem_free=15G,h_vmem=16G,h_fsize=200G -pe local 25

library(devtools)
devtools::install_github("PeteHaitch/bsseq")
library(bsseq)
devtools::session_info()
setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/")
load("BSseq_BA24.rda")
BSseq=BSseq_BA24
rm(BSseq_BA24)
gc()
pData(BSseq)$Sex_NeuN <- paste(pData(BSseq)$Sex, pData(BSseq)$NeuN,sep = "_")
design <- model.matrix(~BSseq$Sex* BSseq$NeuN)
colnames(design) <- gsub("BSseq\\$", "", colnames(design))
contrasts <- diag(rep(1,4))
contrasts <- contrasts[, -1]
rownames(contrasts) <- colnames(design)
fac <- pData(BSseq)$Sex_NeuN

set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq, design = design,
                               contrasts = contrasts,
                               cutoff = c(-4.6 ^ 2, 4.6 ^ 2), fac = fac,
                               nperm = 1000, coef = NULL, maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300, mc.cores = 20)

setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/")
Sex_dmrs_NA=dmrs
Sex_dmrs_HC=dmrs
save(Sex_dmrs_NA, file = "Sex_F-stat_DMRs_NA_nperm1000.rda")
save(Sex_dmrs_HC, file = "Sex_F-stat_DMRs_HC_nperm1000.rda")

bstat=dmrs$bstat

pdf(file = "HC_Sex_F-stats.pdf")
plotManyRegions(BSseq = BSseq,
                BSseqStat = bstat,
                regions = x[12200:12271,],col=pData(BSseq)$Sex_color,
                stat = "stat", addRegions = x, geneTrack=RefSeq_exons3,extend = 5000,
                stat.ylim = c(-10, 10))
dev.off()

###################
#
# NEW WAY FROM PETE

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/BS.fit.small.sorted.somatic.all.rda")

#### HC pos
#qrsh -l cegs,mem_free=20G,h_vmem=21G,h_fsize=200G -pe local 20

BSseq.pos.HC <- BS.fit.small.sorted.somatic.all[, which(BS.fit.small.sorted.somatic.all$Tissue == "HC" & BS.fit.small.sorted.somatic.all$NeuN == "pos")]
design0 <- model.matrix(~0 + BSseq.pos.HC$Sex)
colnames(design0) <- gsub("BSseq.pos.HC\\$Sex", "", colnames(design0))

contrasts0 <- limma::makeContrasts(
  MvsFinHC_pos = male - female,
  levels = design0
)

set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.pos.HC,
                               design = design0,
                               contrasts = contrasts0,
                               cutoff = 4.6,
                               fac = BSseq.pos.HC$Sex,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 20)
save(dmrs, file = "Sex_F-stat_1000perm_HC_pos_dmrs.rda")

#### NA pos
#qrsh -l cegs,mem_free=20G,h_vmem=21G,h_fsize=200G -pe local 20

BSseq.pos.NA <- BS.fit.small.sorted.somatic.all[, which(BS.fit.small.sorted.somatic.all$Tissue == "NA" & BS.fit.small.sorted.somatic.all$NeuN == "pos")]
design_sex0 <- model.matrix(~ 0 + BSseq.pos.NA$Sex)
colnames(design_sex0) <- gsub("BSseq.pos.NA\\$Sex", "",
                              colnames(design_sex0))

contrasts0 <- limma::makeContrasts(
  MvsF_NApos = male - female,
  levels = design_sex0
)
set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.pos.NA,
                               design = design_sex0,
                               contrasts = contrasts0,
                               cutoff = 4.6,
                               fac = BSseq.pos.NA$Sex,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 20)
save(dmrs, file = "Sex_F-stat_1000perm_NA_pos_dmrs.rda")

##### BA9 pos
BSseq.pos.BA9 <- BS.fit.small.sorted.somatic.all[, which(BS.fit.small.sorted.somatic.all$Tissue == "BA9" & BS.fit.small.sorted.somatic.all$NeuN == "pos")]

design_sex0 <- model.matrix(~ 0 + BSseq.pos.BA9$Sex)
colnames(design_sex0) <- gsub("BSseq.pos.BA9\\$Sex", "",
                              colnames(design_sex0))

contrasts0 <- limma::makeContrasts(
  MvsF_BA9pos = male - female,
  levels = design_sex0
)
set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.pos.BA9,
                               design = design_sex0,
                               contrasts = contrasts0,
                               cutoff = 4.6,
                               fac = BSseq.pos.BA9$Sex,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 20)
save(dmrs, file = "Sex_F-stat_1000perm_BA9_pos_dmrs.rda")

##### BA24 pos
BSseq.pos.BA24 <- BS.fit.small.sorted.somatic.all[, which(BS.fit.small.sorted.somatic.all$Tissue == "BA24" & BS.fit.small.sorted.somatic.all$NeuN == "pos")]

design_sex0 <- model.matrix(~ 0 + BSseq.pos.BA24$Sex)
colnames(design_sex0) <- gsub("BSseq.pos.BA24\\$Sex", "",
                              colnames(design_sex0))

contrasts0 <- limma::makeContrasts(
  MvsF_BA24pos = male - female,
  levels = design_sex0
)
set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.pos.BA24,
                               design = design_sex0,
                               contrasts = contrasts0,
                               cutoff = 4.6,
                               fac = BSseq.pos.BA24$Sex,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 20)
save(dmrs, file = "Sex_F-stat_1000perm_BA24_pos_dmrs.rda")

#### Just look at sex diff regardless of tissue (pos only)
# qrsh 80G -pe local 4

BSseq.pos <- BS.fit.small.sorted.somatic.all[, which(BS.fit.small.sorted.somatic.all$NeuN == "pos")]

design_sex0 <- model.matrix(~ 0 + BSseq.pos$Sex)
colnames(design_sex0) <- gsub("BSseq.pos\\$Sex", "",
                              colnames(design_sex0))

contrasts0 <- limma::makeContrasts(
  MvsF_pos = male - female,
  levels = design_sex0
)
set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.pos,
                               design = design_sex0,
                               contrasts = contrasts0,
                               cutoff = 4.6,
                               fac = BSseq.pos$Sex,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 4)

save(dmrs, file = "Sex_F-stat_1000perm_pos_dmrs.rda")

#### Just look at sex diff regardless of tissue (neg only)
# qrsh 80G -pe local 4
BSseq.neg <- BS.fit.small.sorted.somatic.all[, which(BS.fit.small.sorted.somatic.all$NeuN == "neg")]

design_sex0 <- model.matrix(~ 0 + BSseq.neg$Sex)
colnames(design_sex0) <- gsub("BSseq.neg\\$Sex", "",
                              colnames(design_sex0))

contrasts0 <- limma::makeContrasts(
  MvsF_pos = male - female,
  levels = design_sex0
)
set.seed(12345)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.neg,
                               design = design_sex0,
                               contrasts = contrasts0,
                               cutoff = 4.6,
                               fac = BSseq.neg$Sex,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 4)

save(dmrs, file = "Sex_F-stat_1000perm_neg_dmrs.rda")

#-------------------------------------------------------------------------------
# F-stat permutations of NeuN+ samples excluding nucleus accumbens (NA)

pos_nonNA <- intersect(grep("pos", colnames(BSseq)),
                       grep("NA", colnames(BSseq), invert = TRUE))
BSseq.NeuNpos.nonNA <- BSseq[, pos_nonNA]

design <- model.matrix(~BSseq.NeuNpos.nonNA$Tissue)
colnames(design) <- gsub("BSseq.NeuNpos.nonNA\\$Tissue", "", colnames(design))

contrasts <- diag(3)[, -1]

rownames(contrasts) <- colnames(design)
fac <- pData(BSseq.NeuNpos.nonNA)$Tissue

set.seed(666)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.NeuNpos.nonNA,
                               design = design,
                               contrasts = contrasts,
                               cutoff = 4.6 ^ 2,
                               fac = fac,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 6)

save(dmrs,
     file = "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/NeuNpos.nonNA.f-stat_1000perm_dmrs.rda")

pos_nonNA_dmrs <- dmrs$dmrs
sig_pos_nonNA_dmrs <- pos_nonNA_dmrs[which(pos_nonNA_dmrs$fwer <= 50), ]
sig_pos_nonNA_dmrs_gr <- GRanges(sig_pos_nonNA_dmrs) # 208; 347 with fwer 100; 571 with 200 

save(sig_pos_nonNA_dmrs_gr,
     file = "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/All_NeuNpos_nonNA_DMRs_fwer50.rda")


#### Repeat above using 4^2 cutoff rather than 4.6^2
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.small.somatic.all.pos.rda")


pos_nonNA <- grep("NA", colnames(BSseq.pos), invert = TRUE)
BSseq.NeuNpos.nonNA <- BSseq.pos[, pos_nonNA]

design <- model.matrix(~BSseq.NeuNpos.nonNA$Tissue)
colnames(design) <- gsub("BSseq.NeuNpos.nonNA\\$Tissue", "", colnames(design))

contrasts <- diag(3)[, -1]

rownames(contrasts) <- colnames(design)
fac <- pData(BSseq.NeuNpos.nonNA)$Tissue

set.seed(666)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.NeuNpos.nonNA,
                               design = design,
                               contrasts = contrasts,
                               cutoff = 4 ^ 2,
                               fac = fac,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 6)

save(dmrs,
     file = "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/NeuNpos.nonNA.f-stat_cutoff4_1000perm_dmrs.rda")

pos_nonNA_dmrs <- dmrs$dmrs
sig_pos_nonNA_dmrs <- pos_nonNA_dmrs[which(pos_nonNA_dmrs$fwer <= 50), ]
sig_pos_nonNA_dmrs_gr <- GRanges(sig_pos_nonNA_dmrs) # 


sig_pos_nonNA_dmrs <- pos_nonNA_dmrs[which(pos_nonNA_dmrs$fwer <= 500), ]
nrow(sig_pos_nonNA_dmrs)  #1317
sig_pos_nonNA_dmrs_gr <- GRanges(sig_pos_nonNA_dmrs)  
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/RefSeq_exons3.rds")

pdf("nonNA_bottomDMRs.pdf")
plotManyRegions(BSseq = BSseq.NeuNpos.nonNA,
                regions =sig_pos_nonNA_dmrs_gr[1200:1317,],col=pData(BSseq.NeuNpos.nonNA)$Tissue_color,
                stat = "stat", addRegions = sig_pos_nonNA_dmrs_gr, geneTrack=RefSeq_exons3,extend = 5000,
                stat.ylim = c(-10, 10))

dev.off()

pdf("nonNA_topDMRs.pdf")
plotManyRegions(BSseq = BSseq.NeuNpos.nonNA,
                regions =sig_pos_nonNA_dmrs_gr[1:500,],col=pData(BSseq.NeuNpos.nonNA)$Tissue_color,
                stat = "stat", addRegions = sig_pos_nonNA_dmrs_gr, geneTrack=RefSeq_exons3,extend = 5000,
                stat.ylim = c(-10, 10))
dev.off()


# are there any in this groupt that aren't in sig_pos_dmrs?

plot=setdiff(sig_pos_nonNA_dmrs_gr,sig_pos_dmrs_gr)
j=as.data.frame(plot)
j=j[,1]
plot2=sig_pos_nonNA_dmrs_gr[!j]

#####################  repeat with lower cutoff

#### Repeat above using 3^2 cutoff rather than 4^2
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.small.somatic.all.pos.rda")


pos_nonNA <- grep("NA", colnames(BSseq.pos), invert = TRUE)
BSseq.NeuNpos.nonNA <- BSseq.pos[, pos_nonNA]

design <- model.matrix(~BSseq.NeuNpos.nonNA$Tissue)
colnames(design) <- gsub("BSseq.NeuNpos.nonNA\\$Tissue", "", colnames(design))

contrasts <- diag(3)[, -1]

rownames(contrasts) <- colnames(design)
fac <- pData(BSseq.NeuNpos.nonNA)$Tissue

set.seed(666)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.NeuNpos.nonNA,
                               design = design,
                               contrasts = contrasts,
                               cutoff = 3 ^ 2,
                               fac = fac,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 7)

save(dmrs,
     file = "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/NeuNpos.nonNA.f-stat_cutoff3_1000perm_dmrs.rda")

pos_nonNA_dmrs <- dmrs$dmrs
sig_pos_nonNA_dmrs <- pos_nonNA_dmrs[which(pos_nonNA_dmrs$fwer <= 50), ]
sig_pos_nonNA_dmrs_gr <- GRanges(sig_pos_nonNA_dmrs) # 280



#### Repeat above using 3.5^2 cutoff rather than 4^2
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/BS.fit.small.somatic.all.pos.rda")


pos_nonNA <- grep("NA", colnames(BSseq.pos), invert = TRUE)
BSseq.NeuNpos.nonNA <- BSseq.pos[, pos_nonNA]

design <- model.matrix(~BSseq.NeuNpos.nonNA$Tissue)
colnames(design) <- gsub("BSseq.NeuNpos.nonNA\\$Tissue", "", colnames(design))

contrasts <- diag(3)[, -1]

rownames(contrasts) <- colnames(design)
fac <- pData(BSseq.NeuNpos.nonNA)$Tissue

set.seed(666)
dmrs <- bsseq:::fstat.pipeline(BSseq = BSseq.NeuNpos.nonNA,
                               design = design,
                               contrasts = contrasts,
                               cutoff = 3.5 ^ 2,
                               fac = fac,
                               nperm = 1000,
                               coef = NULL,
                               maxGap.sd = 10 ^ 8,
                               maxGap.dmr = 300,
                               mc.cores = 5)

save(dmrs,
     file = "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/F-Stat_Analysis/NeuNpos.nonNA.f-stat_cutoff3.5_1000perm_dmrs.rda")

pos_nonNA_dmrs <- dmrs$dmrs
sig_pos_nonNA_dmrs <- pos_nonNA_dmrs[which(pos_nonNA_dmrs$fwer <= 50), ]
sig_pos_nonNA_dmrs_gr <- GRanges(sig_pos_nonNA_dmrs) # 273





