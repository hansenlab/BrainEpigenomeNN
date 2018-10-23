library(bsseq)
library(minfi)
library(devtools)
devtools::session_info()
setwd("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/")
temp=list.files(pattern="*.fit.small.rda")
temp
load(temp[1])
BS=BS.fit.small
for (i in 2:length(temp)) {
	print(cat(paste0(i,"\n"))) #paste0 uses empty space separator
	load(temp[i])
	tmp=BS.fit.small
	BS=combine(BS,tmp)
}

BS.fit.small.sorted=BS

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Flow_pData.rda")
pData(BS.fit.small.sorted)=pData

save(BS.fit.small.sorted,file="BS.fit.small.sorted.rda")

BS.cov=getCoverage(BS.fit.small.sorted)
colnames(BS.cov)=colnames(BS.fit.small.sorted)


pdf("Sorted.Small.Summary.pdf")
densityPlot(bsseq::getMeth(BS.fit.small.sorted),sampGroups=pData$NeuN)
dev.off()

keepLoci.all <- which(rowSums(BS.cov >= 1) >= 45)

BS.fit.small.sorted.all<-BS.fit.small.sorted[keepLoci.all,]
save(BS.fit.small.sorted.all,file="BS.fit.small.sorted.all.rda")

BS.fit.small.sorted.23.all=chrSelectBSseq(BS.fit.small.sorted.all,seqnames=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
BS.fit.small.sorted.somatic.all=chrSelectBSseq(BS.fit.small.sorted.all,seqnames=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
save(BS.fit.small.sorted.23.all,file="BS.fit.small.sorted.23.all.rda")
save(BS.fit.small.sorted.somatic.all,file="BS.fit.small.sorted.somatic.all.rda")



