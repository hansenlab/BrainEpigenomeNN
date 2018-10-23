
### To make final BS objects to use downstream ###

library(devtools)
install_github("kasperdanielhansen/bsseq")
library(bsseq)
library(minfi)
devtools::session_info()

setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation")
load("BS.fit.small.rda") #includes all chromosomes (M and lambda too)

BS.fit.small.23=chrSelectBSseq(BS.fit.small,seqnames=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
BS.fit.small.somatic=chrSelectBSseq(BS.fit.small.23,seqnames=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))


pData=pData(BS.fit.small.23)
pData$NeuN_color=pData$NeuN
pData$NeuN_color=replace(pData$NeuN_color,pData$NeuN_color=="pos","deepskyblue")
pData$NeuN_color=replace(pData$NeuN_color,pData$NeuN_color=="neg","firebrick3")
pData$Tissue_color=as.character(pData$Tissue)
pData$Tissue_color=replace(pData$Tissue_color,pData$Tissue_color=="NA","chocolate1")
pData$Tissue_color=replace(pData$Tissue_color,pData$Tissue_color=="BA9","deepskyblue")
pData$Tissue_color=replace(pData$Tissue_color,pData$Tissue_color=="BA24","deeppink")
pData$Tissue_color=replace(pData$Tissue_color,pData$Tissue_color=="HC","darkgrey")
pData$Sex_color=pData$Sex
pData$Sex_color=replace(pData$Sex_color,pData$Sex_color=="male","deepskyblue")
pData$Sex_color=replace(pData$Sex_color,pData$Sex_color=="female","firebrick3")
pData$Individual_color=pData$Individual
pData$Individual_color=replace(pData$Individual_color,pData$Individual_color=="5248","aquamarine")
pData$Individual_color=replace(pData$Individual_color,pData$Individual_color=="5284","chocolate1")
pData$Individual_color=replace(pData$Individual_color,pData$Individual_color=="5552","firebrick3")
pData$Individual_color=replace(pData$Individual_color,pData$Individual_color=="5570","purple")
pData$Individual_color=replace(pData$Individual_color,pData$Individual_color=="5569","deepskyblue")
pData$Individual_color=replace(pData$Individual_color,pData$Individual_color=="5628","darkgrey")
pData(BS.fit.small.23)$lty <- ifelse(pData(BS.fit.small.23)$NeuN == "pos", 1, 2)

pData(BS.fit.small.23)=pData

pData(BS.fit.small.somatic)=pData

save(BS.fit.small.somatic,file="BS.fit.small.somatic.rda")
save(BS.fit.small.23,file="BS.fit.small.23.rda")

BS.cov=getCoverage(BS.fit.small.23)
colnames(BS.cov)=colnames(BS.fit.small.23)

keepLoci.all <- which(rowSums(BS.cov[, 1:45] >= 1) >= 45)
length(keepLoci.all)  #23876894

BS.fit.small.all<-BS.fit.small.23[keepLoci.all,]
pData(BS.fit.small.all)=pData
save(BS.fit.small.all,file="BS.fit.small.23.all.rda")

BS.fit.small.somatic.all<-BS.fit.small.somatic[keepLoci.all,]
save(BS.fit.small.somatic.all,file="BS.fit.small.somatic.all.rda")

#####################################################################################
