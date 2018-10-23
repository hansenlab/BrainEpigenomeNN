######  Use same BS.all file for all comparisons! ######


library(devtools)
library(bsseq)
library(minfi)
devtools::session_info()
setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/")
load("BS.fit.small.sorted.somatic.all.rda") 

setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/t-stat_Analysis/")
pData=pData(BS.fit.small.sorted.somatic.all)

names=sampleNames(BS.fit.small.sorted.somatic.all)
pos=names[grepl("pos",names)] #subset just the pos
neg=names[grepl("neg",names)] #subset just the neg

pos.fit.small.all=BS.fit.small.sorted.somatic.all[,pos]
neg.fit.small.all=BS.fit.small.sorted.somatic.all[,neg]
pData(pos.fit.small.all)=pData(BS.fit.small.sorted.somatic.all[,pos])
pData(neg.fit.small.all)=pData(BS.fit.small.sorted.somatic.all[,neg])

NApos=pos[grepl("NA",pos)]
HCpos=pos[grepl("HC",pos)]
NA.HC.pos.tstat<-BSmooth.tstat(pos.fit.small.all,group1=NApos,group2=HCpos, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(NA.HC.pos.tstat,file="NA.HC.pos.tstat.rda")

dmrs0 <- dmrFinder(NA.HC.pos.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs0) 
nrow(dmrs)  
pdf("NA.HC.pos.DMRs.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all[,c(NApos,HCpos)],dmrs[1:100,],extend=10000,addRegions=dmrs,col=pData(pos.fit.small.all[,c(NApos,HCpos)])$Tissue_color,main="NAvsHC.pos")
dev.off()
pdf("NA.HC.pos.DMRs2.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all,dmrs[1:100,],extend=10000,addRegions=dmrs,col=pData(pos.fit.small.all)$Tissue_color,main="NAvsHC.pos")
dev.off()
source("/amber3/feinbergLab/personal/njung/Leukemic_stem_cell/450k/analysis/code/Rafa/cis-genome.R")
write.csv(cbind(dmrs,matchGenes(dmrs,build="hg19")),file="NA.HC.pos.genelist.csv")

# Now do NA BA9
NApos=pos[grepl("NA",pos)]
BA9pos=pos[grepl("BA9",pos)]
NA.BA9.pos.tstat<-BSmooth.tstat(pos.fit.small.all,group1=NApos,group2=BA9pos, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(NA.BA9.pos.tstat,file="NA.BA9.pos.tstat.rda")

dmrs0 <- dmrFinder(NA.BA9.pos.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs1 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs1)
pdf("NA.BA9.pos.DMRs.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all[,c(NApos,BA9pos)],dmrs1[1:100,],extend=10000,addRegions=dmrs1,col=pData(pos.fit.small.all[,c(NApos,BA9pos)])$Tissue_color,main="NAvsBA9.pos")
dev.off()
pdf("NA.BA9.pos.DMRs2.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all,dmrs1[1:100,],extend=10000,addRegions=dmrs1,col=pData(pos.fit.small.all)$Tissue_color,main="NAvsBA9.pos")
dev.off()
write.csv(cbind(dmrs1,matchGenes(dmrs1,build="hg19")),file="NA.BA9.pos.genelist.csv")

#Compare NA and BA24 pos

NApos=pos[grepl("NA",pos)]
BA24pos=pos[grepl("BA24",pos)]

NA.BA24.pos.tstat<-BSmooth.tstat(pos.fit.small.all,group1=NApos,group2=BA24pos, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(NA.BA24.pos.tstat,file="NA.BA24.pos.tstat.rda")

dmrs0 <- dmrFinder(NA.BA24.pos.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs2 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs0) 
nrow(dmrs2)  
pdf("NA.BA24.pos.DMRs.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all[,c(NApos,BA24pos)],dmrs2[1:100,],extend=10000,addRegions=dmrs2,col=pData(pos.fit.small.all[,c(NApos,BA24pos)])$Tissue_color,main="NAvsBA24.pos")
dev.off()
pdf("NA.BA24.pos.DMRs2.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all,dmrs2[1:100,],extend=10000,addRegions=dmrs2,col=pData(pos.fit.small.all)$Tissue_color,main="NAvsBA24.pos")
dev.off()
write.csv(cbind(dmrs2,matchGenes(dmrs2,build="hg19")),file="NA.BA24.pos.genelist.csv")

#HC vs BA9
HCpos=pos[grepl("HC",pos)]
BA9pos=pos[grepl("BA9",pos)]

HC.BA9.pos.tstat<-BSmooth.tstat(pos.fit.small.all,group1=HCpos,group2=BA9pos, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(HC.BA9.pos.tstat,file="HC.BA9.pos.tstat.rda")

dmrs0 <- dmrFinder(HC.BA9.pos.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs3 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs0) 
nrow(dmrs3)  
pdf("HC.BA9.pos.DMRs.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all[,c(BA9pos,HCpos)],dmrs3[1:100,],extend=10000,addRegions=dmrs3,col=pData(pos.fit.small.all[,c(BA9pos,HCpos)])$Tissue_color,main="HCvsBA9.pos")
dev.off()
pdf("HC.BA9.pos.DMRs2.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all,dmrs3[1:100,],extend=10000,addRegions=dmrs3,col=pData(pos.fit.small.all)$Tissue_color,main="HCvsBA9.pos")
dev.off()
write.csv(cbind(dmrs3,matchGenes(dmrs3,build="hg19")),file="HC.BA9.pos.genelist.csv")

#BA24 vs BA9
BA24pos=pos[grepl("BA24",pos)]
BA9pos=pos[grepl("BA9",pos)]

BA24.BA9.pos.tstat<-BSmooth.tstat(pos.fit.small.all,group1=BA24pos,group2=BA9pos, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(BA24.BA9.pos.tstat,file="BA24.BA9.pos.tstat.rda")

dmrs0 <- dmrFinder(BA24.BA9.pos.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs4 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs0) 
nrow(dmrs4)  
pdf("BA24.BA9.pos.DMRs.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all[,c(BA9pos,BA24pos)],dmrs4[1:100,],extend=10000,addRegions=dmrs4,col=pData(pos.fit.small.all[,c(BA9pos,BA24pos)])$Tissue_color,main="BA24vsBA9.pos")
dev.off()
pdf("BA24.BA9.pos.DMRs2.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all,dmrs4[1:100,],extend=10000,addRegions=dmrs4,col=pData(pos.fit.small.all)$Tissue_color,main="BA24vsBA9.pos")
dev.off()
write.csv(cbind(dmrs4,matchGenes(dmrs4,build="hg19")),file="BA24.BA9.pos.genelist.csv")

#Compare HC and BA24 pos


HC.BA24.pos.tstat<-BSmooth.tstat(pos.fit.small.all,group1=HCpos,group2=BA24pos, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(HC.BA24.pos.tstat,file="HC.BA24.pos.tstat.rda")

dmrs0 <- dmrFinder(HC.BA24.pos.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs5 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs0) 
nrow(dmrs5)  
pdf("HC.BA24.pos.DMRs.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all[,c(BA24pos,HCpos)],dmrs5[1:100,],extend=10000,addRegions=dmrs5,col=pData(pos.fit.small.all[,c(BA24pos,HCpos)])$Tissue_color,main="HCvsBA24.pos")
dev.off()
pdf("HC.BA24.pos.DMRs2.pdf",width=10,height=5)
plotManyRegions(pos.fit.small.all,dmrs5[1:100,],extend=10000,addRegions=dmrs5,col=pData(pos.fit.small.all)$Tissue_color,main="HCvsBA24.pos")
dev.off()
write.csv(cbind(dmrs5,matchGenes(dmrs5,build="hg19")),file="HC.BA24.pos.genelist.csv")

#Plot all corrected tstats
tstat1 <- getStats(HC.BA24.pos.tstat)[, "tstat.corrected"]
tstat2 <- getStats(BA24.BA9.pos.tstat)[, "tstat.corrected"]
tstat3 <- getStats(HC.BA9.pos.tstat)[, "tstat.corrected"]
tstat4 <- getStats(NA.BA24.pos.tstat)[, "tstat.corrected"]
tstat5 <- getStats(NA.BA9.pos.tstat)[, "tstat.corrected"]
tstat6 <- getStats(NA.HC.pos.tstat)[, "tstat.corrected"]

pdf("All NeuN Pos Tissue Corrected t-stats.pdf")
    plot(density(tstat1), xlim = c(-10,10), ylim=c(0,6),col = "deepskyblue", main = "")
    lines(density(tstat2), col = "darkgrey")
    lines(density(tstat3), col = "chocolate1")
    lines(density(tstat4), col = "deeppink") 
    lines(density(tstat5), col = "green")
    lines(density(tstat6), col = "black") 
    legend("topleft", legend = c("HCvsBA24","BA24vsBA9", "HCvsBA9","NAvsBA24","NAvsBA9","NAvsHC"), lty = c(1,1),
               col = c("deepskyblue", "darkgrey","chocolate1","deeppink","green","black"))
dev.off()

#Plot all uncorrected tstats
tstat1 <- getStats(HC.BA24.pos.tstat)[, "tstat"]
tstat2 <- getStats(BA24.BA9.pos.tstat)[, "tstat"]
tstat3 <- getStats(HC.BA9.pos.tstat)[, "tstat"]
tstat4 <- getStats(NA.BA24.pos.tstat)[, "tstat"]
tstat5 <- getStats(NA.BA9.pos.tstat)[, "tstat"]
tstat6 <- getStats(NA.HC.pos.tstat)[, "tstat"]

pdf("All NeuN Pos Tissue t-stats (uncorrected).pdf")
    plot(density(tstat1), xlim = c(-10,10), ylim=c(0,6), col = "deepskyblue", main = "")
    lines(density(tstat2), col = "darkgrey")
    lines(density(tstat3), col = "chocolate1")
    lines(density(tstat4), col = "deeppink") 
    lines(density(tstat5), col = "green")
    lines(density(tstat6), col = "black") 
    legend("topleft", legend = c("HCvsBA24","BA24vsBA9", "HCvsBA9","NAvsBA24","NAvsBA9","NAvsHC"), lty = c(1,1),
               col = c("deepskyblue", "darkgrey","chocolate1","deeppink","green","black"))
dev.off()


#NOW DO NEGATIVES
rm(pos.fit.small.all)

NAneg=neg[grepl("NA",neg)]
HCneg=neg[grepl("HC",neg)]

NA.HC.neg.tstat<-BSmooth.tstat(neg.fit.small.all,group1=NAneg,group2=HCneg, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(NA.HC.neg.tstat,file="NA.HC.neg.tstat.rda")

dmrs0 <- dmrFinder(NA.HC.neg.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs)
pdf("NA.HC.neg.DMRs.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all[,c(NAneg,HCneg)],dmrs[1:100,],extend=10000,addRegions=dmrs,col=pData(neg.fit.small.all[,c(NAneg,HCneg)])$Tissue_color,main="NAvsHC.neg")
dev.off()
pdf("NA.HC.neg.DMRs2.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all,dmrs[1:100,],extend=10000,addRegions=dmrs,col=pData(neg.fit.small.all)$Tissue_color,main="NAvsHC.neg")
dev.off()
write.csv(cbind(dmrs,matchGenes(dmrs,build="hg19")),file="NAvsHC.neg.genelist.csv")

#compare BA9 and NA neg
BA9neg=neg[grepl("BA9",neg)]

NA.BA9.neg.tstat<-BSmooth.tstat(neg.fit.small.all,group1=NAneg,group2=BA9neg, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(NA.BA9.neg.tstat,file="NA.BA9.neg.tstat.rda")

dmrs0 <- dmrFinder(NA.BA9.neg.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs1 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs1)
pdf("NA.BA9.neg.DMRs.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all[,c(NAneg,BA9neg)],dmrs1[1:100,],extend=10000,addRegions=dmrs1,col=pData(neg.fit.small.all[,c(NAneg,BA9neg)])$Tissue_color,main="NAvsBA9.neg")
dev.off()
pdf("NA.BA9.neg.DMRs2.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all,dmrs1[1:100,],extend=10000,addRegions=dmrs1,col=pData(neg.fit.small.all)$Tissue_color,main="NAvsBA9.neg")
dev.off()
write.csv(cbind(dmrs1,matchGenes(dmrs1,build="hg19")),file="NAvsBA9.neg.genelist.csv")

#NA vs BA24 neg
BA24neg=neg[grepl("BA24",neg)]

NA.BA24.neg.tstat<-BSmooth.tstat(neg.fit.small.all,group1=NAneg,group2=BA24neg, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(NA.BA24.neg.tstat,file="NA.BA24.neg.tstat.rda")

dmrs0 <- dmrFinder(NA.BA24.neg.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs2 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs2)
pdf("NA.BA24.neg.DMRs.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all[,c(NAneg,BA24neg)],dmrs2[1:100,],extend=10000,addRegions=dmrs2,col=pData(neg.fit.small.all[,c(NAneg,BA24neg)])$Tissue_color,main="NAvsBA24.neg")
dev.off()
pdf("NA.BA24.neg.DMRs2.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all,dmrs2[1:100,],extend=10000,addRegions=dmrs2,col=pData(neg.fit.small.all)$Tissue_color,main="NAvsBA24.neg")
dev.off()
write.csv(cbind(dmrs2,matchGenes(dmrs2,build="hg19")),file="NAvsBA24.neg.genelist.csv")

#BA9 vs BA24
BA9.BA24.neg.tstat<-BSmooth.tstat(neg.fit.small.all,group1=BA9neg,group2=BA24neg, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(BA9.BA24.neg.tstat,file="BA9.BA24.neg.tstat.rda")

dmrs0 <- dmrFinder(BA9.BA24.neg.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs3 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs3)
pdf("BA9.BA24.neg.DMRs.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all[,c(BA9neg,BA24neg)],dmrs3[1:100,],extend=10000,addRegions=dmrs3,col=pData(neg.fit.small.all[,c(BA9neg,BA24neg)])$Tissue_color,main="BA9vsBA24.neg")
dev.off()
pdf("BA9.BA24.neg.DMRs2.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all,dmrs3[1:100,],extend=10000,addRegions=dmrs3,col=pData(neg.fit.small.all)$Tissue_color,main="BA9vsBA24.neg")
dev.off()
source("/amber3/feinbergLab/personal/njung/Leukemic_stem_cell/450k/analysis/code/Rafa/cis-genome.R")
write.csv(cbind(dmrs3,matchGenes(dmrs3,build="hg19")),file="BA9vsBA24.neg.genelist.csv")

#compare BA9 and HC neg
BA9.HC.neg.tstat<-BSmooth.tstat(neg.fit.small.all,group1=BA9neg,group2=HCneg, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(BA9.HC.neg.tstat,file="BA9.HC.neg.tstat.rda")

dmrs0 <- dmrFinder(BA9.HC.neg.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs4 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs0) #256767
nrow(dmrs4)  #146228
pdf("BA9.HC.neg.DMRs.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all[,c(BA9neg,HCneg)],dmrs4[1:100,],extend=10000,addRegions=dmrs4,col=pData(neg.fit.small.all[,c(BA9neg,HCneg)])$Tissue_color,main="BA9vsHC.neg")
dev.off()
pdf("BA9.HC.neg.DMRs2.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all,dmrs4[1:100,],extend=10000,addRegions=dmrs4,col=pData(neg.fit.small.all)$Tissue_color,main="BA9vsHC.neg")
dev.off()
write.csv(cbind(dmrs4,matchGenes(dmrs4,build="hg19")),file="BA9vsHC.neg.genelist.csv")


#compare BA24 and HC neg
BA24.HC.neg.tstat<-BSmooth.tstat(neg.fit.small.all,group1=BA24neg,group2=HCneg, estimate.var="same", local.correct=TRUE, verbose=TRUE,maxGap=300)
save(BA24.HC.neg.tstat,file="BA24.HC.neg.tstat.rda")

dmrs0 <- dmrFinder(BA24.HC.neg.tstat, cutoff=c(-4.6,4.6), stat="tstat.corrected") 
dmrs5 <- dmrs0[which(dmrs0$n>=5 & abs(dmrs0$meanDiff)>0.1),]
nrow(dmrs0) #256767
nrow(dmrs5)  #146228
pdf("BA24.HC.neg.DMRs.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all[,c(BA24neg,HCneg)],dmrs5[1:100,],extend=10000,addRegions=dmrs5,col=pData(neg.fit.small.all[,c(BA24neg,HCneg)])$Tissue_color,main="BA24vsHC.neg")
dev.off()
pdf("BA24.HC.neg.DMRs2.pdf",width=10,height=5)
plotManyRegions(neg.fit.small.all,dmrs5[1:100,],extend=10000,addRegions=dmrs5,col=pData(neg.fit.small.all)$Tissue_color,main="BA24vsHC.neg")
dev.off()
write.csv(cbind(dmrs5,matchGenes(dmrs5,build="hg19")),file="BA24vsHC.neg.genelist.csv")


#Plot all corrected tstats
tstat1 <- getStats(BA24.HC.neg.tstat)[, "tstat.corrected"]
tstat2 <- getStats(BA9.BA24.neg.tstat)[, "tstat.corrected"]
tstat3 <- getStats(BA9.HC.neg.tstat)[, "tstat.corrected"]
tstat4 <- getStats(NA.BA24.neg.tstat)[, "tstat.corrected"]
tstat5 <- getStats(NA.BA9.neg.tstat)[, "tstat.corrected"]
tstat6 <- getStats(NA.HC.neg.tstat)[, "tstat.corrected"]

pdf("All NeuN Neg Tissue Corrected t-stats.pdf")
    plot(density(tstat1), xlim = c(-10,10),ylim=c(0,6), col = "deepskyblue", main = "")
    lines(density(tstat2), col = "darkgrey")
    lines(density(tstat3), col = "chocolate1")
    lines(density(tstat4), col = "deeppink") 
    lines(density(tstat5), col = "green")
    lines(density(tstat6), col = "black") 
    legend("topleft", legend = c("HCvsBA24","BA24vsBA9", "HCvsBA9","NAvsBA24","NAvsBA9","NAvsHC"), lty = c(1,1),
               col = c("deepskyblue", "darkgrey","chocolate1","deeppink","green","black"))
dev.off()

#Plot uncorrected tstats
tstat1 <- getStats(BA24.HC.neg.tstat)[, "tstat"]
tstat2 <- getStats(BA9.BA24.neg.tstat)[, "tstat"]
tstat3 <- getStats(BA9.HC.neg.tstat)[, "tstat"]
tstat4 <- getStats(NA.BA24.neg.tstat)[, "tstat"]
tstat5 <- getStats(NA.BA9.neg.tstat)[, "tstat"]
tstat6 <- getStats(NA.HC.neg.tstat)[, "tstat"]

pdf("All NeuN Neg Tissue t-stats.pdf")
    plot(density(tstat1), xlim = c(-10,10),ylim=c(0,6), col = "deepskyblue", main = "")
    lines(density(tstat2), col = "darkgrey")
    lines(density(tstat3), col = "chocolate1")
    lines(density(tstat4), col = "deeppink") 
    lines(density(tstat5), col = "green")
    lines(density(tstat6), col = "black") 
    legend("topleft", legend = c("HCvsBA24","BA24vsBA9", "HCvsBA9","NAvsBA24","NAvsBA9","NAvsHC"), lty = c(1,1),
               col = c("deepskyblue", "darkgrey","chocolate1","deeppink","green","black"))
dev.off()

