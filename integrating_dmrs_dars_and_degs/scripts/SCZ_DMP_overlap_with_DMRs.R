library(GenomicRanges)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
x <- read.csv("https://images.nature.com/original/nature-assets/neuro/journal/v19/n1/extref/nn.4181-S12.csv")
x_gr <- makeGRangesFromDataFrame(Locations[x$X, ], start.field = "pos", seqnames.field = "chr", end.field = "pos")

load("FlowSortingProject/integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")


table(overlapsAny(x_gr, dmrs_pos, ignore.strand = TRUE))
table(overlapsAny(x_gr, Annotated_POSvNEG_DMRs_gr))
