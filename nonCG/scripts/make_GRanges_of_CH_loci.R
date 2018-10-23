# Construct GRanges with CA and CT loci in hg19
# Peter Hickey
# 2017-10-12

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

extdir <- "../extdata"

pos_CA <- bsapply(new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19,
                      FUN = function(s) matchPattern("CA", s)))
pos_CA_gr <- resize(GRanges(seqnames = Rle(names(pos_CA), lengths(pos_CA)),
                            ranges = unlist(as(pos_CA, "IRangesList"),
                                            use.names = FALSE),
                            strand = "+",
                            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
                    width = 1,
                    fix = "start")
saveRDS(pos_CA_gr, file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                          "pos_CA.loci.GRanges.rds"))

neg_CA <- bsapply(new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19,
                      FUN = function(s) matchPattern("TG", s)))
neg_CA_gr <- resize(GRanges(seqnames = Rle(names(neg_CA), lengths(neg_CA)),
                            ranges = unlist(as(neg_CA, "IRangesList"),
                                            use.names = FALSE),
                            strand = "-",
                            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
                    width = 1,
                    fix = "start")
saveRDS(neg_CA_gr, file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                             "neg_CA.loci.GRanges.rds"))

pos_CT <- bsapply(new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19,
                      FUN = function(s) matchPattern("CT", s)))
pos_CT_gr <- resize(GRanges(seqnames = Rle(names(pos_CT), lengths(pos_CT)),
                            ranges = unlist(as(pos_CT, "IRangesList"),
                                            use.names = FALSE),
                            strand = "+",
                            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
                    width = 1,
                    fix = "start")
saveRDS(pos_CT_gr, file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                             "pos_CT.loci.GRanges.rds"))

neg_CT <- bsapply(new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19,
                      FUN = function(s) matchPattern("AG", s)))
neg_CT_gr <- resize(GRanges(seqnames = Rle(names(neg_CT), lengths(neg_CT)),
                            ranges = unlist(as(neg_CT, "IRangesList"),
                                            use.names = FALSE),
                            strand = "-",
                            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
                    width = 1,
                    fix = "start")
saveRDS(neg_CT_gr, file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                             "neg_CT.loci.GRanges.rds"))

