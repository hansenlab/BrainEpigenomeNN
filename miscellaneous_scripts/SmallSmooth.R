library(bsseq)
library(devtools)
devtools::session_info()
args=(commandArgs(TRUE))
print (args[1])
print (args[2])
file=args[1]
name=strsplit(file,'[.]')
name=sapply(name,"[[",1)
BS = read.bismark(file,name,fileType="cytosineReport",rmZeroCov = FALSE, verbose = TRUE)
BS=sort(sortSeqlevels(BS))
BS.fit.small<-BSmooth(BS,verbose=TRUE,ns=20,h=1000)
save(BS.fit.small, file = paste(name,".BS.fit.small.rda",sep=""),compress=T)

