### EBV.csv from Hansen et al 2014- SupTable2_SmallDMRs - Transformation tab

EBV=read.csv("/Users/Lindsay/Desktop/EBV.csv",header=T)
EBV_gr=GRanges(EBV)
load("/Volumes/Transcend/assays-and-features.rda")
source("/Users/Lindsay/Desktop/GTExScripts/FlowSortingProject/integrating-dmrs-dars-and-degs/scripts/functions.R")
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(stringr)
x=readRDS("/Users/Lindsay/Desktop/GTExScripts/FlowSortingProject/integrating-dmrs-dars-and-degs/objects/CG_DMRs_OR_df.rds")


list_of_CG_DMRs <- c(list(
        "CG-DMRs (POS)" = dmrs_pos,
        "CG-DMRs (POSvsNEG)" = Annotated_POSvNEG_DMRs_gr, "EBV-DMRs"=EBV_gr))
non_degs <- unflattened_features$genes[non_deg_names]
deg_promoters <- promoters[match(deg_names, unlist(promoters$GENEID))]
non_deg_promoters <- promoters[match(non_deg_names, unlist(promoters$GENEID))]


list_of_OCRs <- list(
        "union" = ocrs_overall)

features <- c(
        cgi_features,
        gencode_features[["pc_transcripts"]][c("promoter", "five_utr", "exonic",
                                               "intronic", "three_utr",
                                               "intergenic")],
        unlinked_enhancers[c("H3K27ac", "FANTOM5")],
        list("OCR (union)" = list_of_OCRs[["union"]]),
        list_of_CG_DMRs[c("CG-DMRs (POS)", "CG-DMRs (POSvsNEG)")])


CG_DMRs_OR_df <- bind_rows(
        lapply(names(list_of_CG_DMRs), function(name) {
                CG_DMRs <- list_of_CG_DMRs[[name]]
                bind_rows(FT(x = subsetByOverlaps(cpgs, CG_DMRs),
                             notx = subsetByOverlaps(cpgs, CG_DMRs, invert = TRUE),
                             features = features,
                             db = name),
                          FT(x = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs),
                                                  c(degs, non_degs, ignore.mcols = TRUE)),
                             notx = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs,
                                                                      invert = TRUE),
                                                     c(degs, non_degs, ignore.mcols = TRUE)),
                             features = list("DEGs" = degs),
                             db = name),
                          FT(x = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs),
                                                  c(deg_promoters, non_deg_promoters,
                                                    ignore.mcols = TRUE)),
                             notx = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs,
                                                                      invert = TRUE),
                                                     c(deg_promoters, non_deg_promoters,
                                                       ignore.mcols = TRUE)),
                             features = list("DEG promoters" = deg_promoters),
                             db = name))
        }))



mCG_OR_df <- bind_rows(CG_DMRs_OR_df) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


hmcol <- colorRampPalette(c("blue", "white", "red"))(100)

mCG_OR_matrix_all <- mCG_OR_df %>%
        dplyr::select(-db, -`CG-DMRs (POS)`, -`CG-DMRs (POSvsNEG)`) %>%
        as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]


heatmap.2(t(mCG_OR_matrix_all),
          Rowv = NULL,
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_OR_matrix_all),
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_OR_matrix_all),
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

