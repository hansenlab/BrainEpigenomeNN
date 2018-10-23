# Enrichment plots for paper
# 2017-02-20
# Peter Hickey

library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(stringr)

load("../../FlowSortingProject/integrating-dmrs-daps-and-degs/objects/assays-and-features.rda")
source("../../FlowSortingProject/integrating-dmrs-daps-and-degs/scripts/functions.R")
cols <- c("E068" = "chocolate1", "E069" = "deeppink", "E071" = "darkgrey",
          "E073" = "deepskyblue")
load("../../FlowSortingProject/integrating-dmrs-daps-and-degs/objects/assays-and-features.rda")

#===============================================================================
# GENCODE and enhancer enrichment
#

#-------------------------------------------------------------------------------
# DMR-CpGs vs. non-DMR-Cpgs: exonic (PC, lncRNA), 5' UTR (PC),
#                           intronic (PC, lncRNA), promoter (PC, lncRNA),
#                           3' UTR (PC), intergenic (union), and enhancers


names(unlinked_enhancers) <- c("Vermunt et al.", "FANTOM5", "Overlap")

or_dmrs_pos_cpgs_gencode_and_enhancers <- bind_rows(
  FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs, gencode_features$union, "union"),
  FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs, gencode_features$pc_transcripts, "PC"),
  FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs, gencode_features$lnc_transcripts,
     "lncRNA"),
  FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs,
     unlinked_enhancers,
     rep("Enhancers", 3))) %>%
     as_data_frame()

g <- or_dmrs_pos_cpgs_gencode_and_enhancers %>%
  filter((db == "PC" & feature %in% c("exonic", "five_utr", "intronic",
                                      "promoter", "three_utr")) |
           (db == "lncRNA" & feature %in% c("exonic", "intronic", "promoter")) |
           (db == "union" & feature == "intergenic") |
           (db == "Enhancers")) %>%
  mutate(feature = case_when(.$feature == "five_utr" ~ "5' UTR",
                             .$feature == "three_utr" ~ "3' UTR",
                             grepl("FANTOM", .$feature) ~ .$feature,
                             .$feature == "Vermunt et al." ~ .$feature,
                             TRUE ~ str_to_title(.$feature)),
         feature = factor(feature,
                          levels = c("Intergenic", "Promoter", "5' UTR",
                                     "Exonic", "Intronic", "3' UTR",
                                     "Vermunt et al.", "FANTOM5",
                                     "Overlap"),
                          ordered = TRUE)) %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("DMR-CpGs vs. non-DMR-CpGs") +
  ylab("log2(OR) with 95% CI") +
  scale_color_discrete(guide = guide_legend(title = "Database")) +
  xlab("Feature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("DMR-CpGs.GENCODE_and_enhancer_enrichment.pdf", plot = g, width = 5,
       height = 5)

#-------------------------------------------------------------------------------
# ATAC peaks vs. rest of genome: exonic (PC, lncRNA), 5' UTR (PC),
#                               intronic (PC, lncRNA), promoter (PC, lncRNA),
#                               3' UTR (PC), intergenic (union), and enhancers

or_bp_peaks_gencode_and_enhancers <- bind_rows(
  FT2(atac_peaks, gencode_features$union, "union", sl),
  FT2(atac_peaks, gencode_features$pc_transcripts, "PC", sl),
  FT2(atac_peaks, gencode_features$lnc_transcripts, "lncRNA", sl),
  FT2(atac_peaks, unlinked_enhancers, rep("Enhancers", 3), sl)) %>%
  as_data_frame()

g <- or_bp_peaks_gencode_and_enhancers %>%
  filter((db == "PC" & feature %in% c("exonic", "five_utr", "intronic",
                                      "promoter", "three_utr")) |
           (db == "lncRNA" & feature %in% c("exonic", "intronic", "promoter")) |
           (db == "union" & feature == "intergenic") |
           (db == "Enhancers")) %>%
  mutate(feature = case_when(.$feature == "five_utr" ~ "5' UTR",
                             .$feature == "three_utr" ~ "3' UTR",
                             grepl("FANTOM", .$feature) ~ .$feature,
                             .$feature == "Vermunt et al." ~ .$feature,
                             TRUE ~ str_to_title(.$feature)),
         feature = factor(feature,
                          levels = c("Intergenic", "Promoter", "5' UTR",
                                     "Exonic", "Intronic", "3' UTR",
                                     "Vermunt et al.", "FANTOM5",
                                     "Overlap"),
                          ordered = TRUE)) %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Base-level: ATAC-peaks vs. rest of genome") +
  ylab("log2(OR) with 95% CI") +
  scale_color_discrete(guide = guide_legend(title = "Database")) +
  xlab("Feature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("ATAC-peaks.GENCODE_and_enhancers_enrichment.pdf", plot = g, width = 5,
       height = 5)

#-------------------------------------------------------------------------------
# ATAC peaks vs. rest of genome (DEGs): exonic (PC, lncRNA), 5' UTR (PC),
#                                       intronic (PC, lncRNA),
#                                       promoter (PC, lncRNA), 3' UTR (PC),
#                                       intergenic (union), and enhancers

or_bp_peaks_vs_rest_of_genome_gencode_degs <- bind_rows(
  FT2(atac_peaks, deg_flattened_features_union, "DEGs-union", sl),
  FT2(atac_peaks, deg_flattened_features_pc, "DEGs-PC", sl),
  FT2(atac_peaks, deg_flattened_features_lnc, "DEGs-lncRNA", sl)) %>%
  as_data_frame()

g <- or_bp_peaks_vs_rest_of_genome_gencode_degs %>%
  filter(feature != "genic") %>%
  mutate(feature = case_when(.$feature == "five_utr" ~ "5' UTR",
                             .$feature == "three_utr" ~ "3' UTR",
                             grepl("FANTOM", .$feature) ~ .$feature,
                             .$feature == "Vermunt et al." ~ .$feature,
                             TRUE ~ str_to_title(.$feature)),
         feature = factor(feature,
                          levels = c("Intergenic", "Promoter", "5' UTR",
                                     "Exonic", "Intronic", "3' UTR",
                                     "Vermunt et al.", "FANTOM5",
                                     "Overlap"),
                          ordered = TRUE)) %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Base-level: ATAC-peaks vs. rest of genome (DEGs)") +
  ylab("log2(OR) with 95% CI") +
  scale_color_discrete(guide = guide_legend(title = "Database")) +
  xlab("Feature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("ATAC-peaks.GENCODE_and_enhancers_enrichment.DEGs.pdf",
       plot = g,
       width = 5,
       height = 5)

#-------------------------------------------------------------------------------
# DAPs vs. rest of genome: exonic (PC, lncRNA), 5' UTR (PC),
#                          intronic (PC, lncRNA), promoter (PC, lncRNA),
#                          3' UTR (PC), intergenic (union), and enhancers
#

or_bp_daps_vs_rest_of_genome_gencode_and_enhancers <- bind_rows(
  FT2(daps, gencode_features$union, "union", sl),
  FT2(daps, gencode_features$pc_transcripts, "PC", sl),
  FT2(daps, gencode_features$lnc_transcripts, "lncRNA", sl),
  FT2(daps, unlinked_enhancers, rep("Enhancers", 3), sl)) %>%
  as_data_frame()

g <- or_bp_daps_vs_rest_of_genome_gencode_and_enhancers %>%
  filter((db == "PC" & feature %in% c("exonic", "five_utr", "intronic",
                                      "promoter", "three_utr")) |
           (db == "lncRNA" & feature %in% c("exonic", "intronic", "promoter")) |
           (db == "union" & feature == "intergenic") |
           (db == "Enhancers")) %>%
  mutate(feature = case_when(.$feature == "five_utr" ~ "5' UTR",
                             .$feature == "three_utr" ~ "3' UTR",
                             grepl("FANTOM", .$feature) ~ .$feature,
                             .$feature == "Vermunt et al." ~ .$feature,
                             TRUE ~ str_to_title(.$feature)),
         feature = factor(feature,
                          levels = c("Intergenic", "Promoter", "5' UTR",
                                     "Exonic", "Intronic", "3' UTR",
                                     "Vermunt et al.", "FANTOM5",
                                     "Overlap"),
                          ordered = TRUE)) %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Base-level: DAPs vs. rest of genome") +
  ylab("log2(OR) with 95% CI") +
  scale_color_discrete(guide = guide_legend(title = "Database")) +
  xlab("Feature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("DAPs_vs_rest_of_genome.GENCODE_and_enhancers_enrichment.pdf",
       plot = g,
       width = 5,
       height = 5)

#-------------------------------------------------------------------------------
# DAPs vs. rest of genome (DEGs): exonic (PC, lncRNA), 5' UTR (PC),
#                                 intronic (PC, lncRNA), promoter (PC, lncRNA),
#                                 3' UTR (PC), intergenic (union), and enhancers

or_bp_daps_vs_rest_of_genome_gencode_degs <- bind_rows(
  FT2(daps, deg_flattened_features_union, "DEGs-union", sl),
  FT2(daps, deg_flattened_features_pc, "DEGs-PC", sl),
  FT2(daps, deg_flattened_features_lnc, "DEGs-lncRNA", sl)) %>%
  as_data_frame()

g <- or_bp_daps_vs_rest_of_genome_gencode_degs %>%
  filter(feature != "genic") %>%
  mutate(feature = case_when(.$feature == "five_utr" ~ "5' UTR",
                             .$feature == "three_utr" ~ "3' UTR",
                             TRUE ~ str_to_title(.$feature)),
         feature = factor(feature,
                          levels = c("Intergenic", "Promoter", "5' UTR",
                                     "Exonic", "Intronic", "3' UTR"),
                          ordered = TRUE)) %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Base-level: DAPs vs. rest of genome (DEGs)") +
  ylab("log2(OR) with 95% CI") +
  scale_color_discrete(guide = guide_legend(title = "Database")) +
  xlab("Feature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("DAPs_vs_rest_of_genome.GENCODE_and_enhancers_enrichment.DEGs.pdf",
       plot = g,
       width = 5,
       height = 5)

#-------------------------------------------------------------------------------
# DAPs vs. non-DAPs: exonic (PC, lncRNA), 5' UTR (PC), intronic (PC, lncRNA),
#                   promoter (PC, lncRNA), 3' UTR (PC), intergenic (union),
#                   and enhancers
#

or_bp_daps_vs_non_daps_gencode_and_enhancers <- bind_rows(
  FT3(daps, non_daps, gencode_features$union, "union"),
  FT3(daps, non_daps, gencode_features$pc_transcripts, "PC"),
  FT3(daps, non_daps, gencode_features$lnc_transcripts, "lncRNA"),
  FT3(daps, non_daps, unlinked_enhancers, rep("Enhancers", 3))) %>%
  as_data_frame()

g <- or_bp_daps_vs_non_daps_gencode_and_enhancers %>%
  filter((db == "PC" & feature %in% c("exonic", "five_utr", "intronic",
                                      "promoter", "three_utr")) |
           (db == "lncRNA" & feature %in% c("exonic", "intronic", "promoter")) |
           (db == "union" & feature == "intergenic") |
           (db == "Enhancers")) %>%
  mutate(feature = case_when(.$feature == "five_utr" ~ "5' UTR",
                             .$feature == "three_utr" ~ "3' UTR",
                             grepl("FANTOM", .$feature) ~ .$feature,
                             .$feature == "Vermunt et al." ~ .$feature,
                             TRUE ~ str_to_title(.$feature)),
         feature = factor(feature,
                          levels = c("Intergenic", "Promoter", "5' UTR",
                                     "Exonic", "Intronic", "3' UTR",
                                     "Vermunt et al.", "FANTOM5",
                                     "Overlap"),
                          ordered = TRUE)) %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Base-level: DAPs vs. non-DAPs") +
  ylab("log2(OR) with 95% CI") +
  scale_color_discrete(guide = guide_legend(title = "Database")) +
  xlab("Feature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("DAPs_vs_non-DAPs.GENCODE_and_enhancers_enrichment.pdf",
       plot = g,
       width = 5,
       height = 5)

#===============================================================================
# CGI enrichment
#

#-------------------------------------------------------------------------------
# DMR-CpGs vs. non-DMR-CpGs: CGIs, shores, shelves, open sea
#

or_dmrs_pos_cpgs_cgi <- FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs, cgi_features,
                           "UCSC")
g <- or_dmrs_pos_cpgs_cgi %>%
  mutate(feature = gsub("CGI", "CpG islands", gsub("OpenSea", "Open sea",
                                                  feature)),
    feature = factor(feature,
                          levels = c("CpG islands", "Shores", "Shelves",
                                     "Open sea"),
                     ordered = TRUE)) %>%
  ggplot(aes(x = feature, y = log2(estimate))) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("DMR-CpGs vs. non-DMR-CpGs") +
  ylab("log2(OR) with 95% CI") +
  scale_color_discrete(guide = guide_legend(title = "Tx db")) +
  xlab("Feature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("DMR-CpGs.CGI-enrichment.pdf", plot = g, width = 5, height = 5)

#===============================================================================
# chromHMM enrichment
#

#-------------------------------------------------------------------------------
# DMR-CpGs vs. non-DMR-CpGs: chromHMM 15-state model for 4 brain regions
#

or_dmrs_pos_cpgs_AH46921 <- FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs, AH46921,
                               "E068")
or_dmrs_pos_cpgs_AH46922 <- FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs, AH46922,
                               "E069")
or_dmrs_pos_cpgs_AH46924 <- FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs, AH46924,
                               "E071")
or_dmrs_pos_cpgs_AH46926 <- FT(dmrs_pos_cpgs, non_dmrs_pos_cpgs, AH46926,
                               "E073")

# NOTE: Order chromHMM states as in
#       http://egg2.wustl.edu/roadmap/figures/mainFigs/Figure_4.jpg
or_dmrs_pos_cpgs_chromhmm <- bind_rows(or_dmrs_pos_cpgs_AH46921,
                                       or_dmrs_pos_cpgs_AH46922,
                                       or_dmrs_pos_cpgs_AH46924,
                                       or_dmrs_pos_cpgs_AH46926) %>%
  mutate(feature = factor(feature,
                          levels = c("Active TSS", "Flanking Active TSS",
                                     "Transcr. at gene 5' and 3'",
                                     "Strong transcription",
                                     "Weak transcription", "Genic enhancers",
                                     "Enhancers", "ZNF genes & repeats",
                                     "Heterochromatin", "Bivalent/Poised TSS",
                                     "Flanking Bivalent TSS/Enh",
                                     "Bivalent Enhancer", "Repressed PolyComb",
                                     "Weak Repressed PolyComb",
                                     "Quiescent/Low"),
                          ordered = TRUE))

g <- or_dmrs_pos_cpgs_chromhmm %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("DMR-CpGs vs. non-DMR-CpGs") +
  ylab("log2(OR) with 95% CI") +
  scale_colour_manual(values = cols)

ggsave("DMR-CpGs.chromHMM-enrichment.pdf", plot = g, width = 5, height = 5)


#-------------------------------------------------------------------------------
# ATAC peaks vs. rest of genome: chromHMM 15-state model for 4 brain regions
#

or_bp_peaks_AH46921 <- FT2(atac_peaks, AH46921, "E068", sl)
or_bp_peaks_AH46922 <- FT2(atac_peaks, AH46922, "E069", sl)
or_bp_peaks_AH46924 <- FT2(atac_peaks, AH46924, "E071", sl)
or_bp_peaks_AH46926 <- FT2(atac_peaks, AH46926, "E073", sl)

or_bp_peaks_chromHMM <- bind_rows(or_bp_peaks_AH46921,
                                  or_bp_peaks_AH46922,
                                  or_bp_peaks_AH46924,
                                  or_bp_peaks_AH46926) %>%
  mutate(feature = factor(feature,
                          levels = c("Active TSS", "Flanking Active TSS",
                                     "Transcr. at gene 5' and 3'",
                                     "Strong transcription",
                                     "Weak transcription", "Genic enhancers",
                                     "Enhancers", "ZNF genes & repeats",
                                     "Heterochromatin", "Bivalent/Poised TSS",
                                     "Flanking Bivalent TSS/Enh",
                                     "Bivalent Enhancer", "Repressed PolyComb",
                                     "Weak Repressed PolyComb",
                                     "Quiescent/Low"),
                          ordered = TRUE))

g <- or_bp_peaks_chromHMM %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Base-level: ATAC-peaks vs. rest of genome") +
  ylab("log2(OR) with 95% CI") +
  scale_colour_manual(values = cols)

ggsave("ATAC-peaks.chromHMM-enrichment.pdf", plot = g, width = 5, height = 5)


#-------------------------------------------------------------------------------
# DAPs vs. non-DAPs (ignoring direction): chromHMM 15-state model for 4 brain
#                                         regions
#

or_bp_daps_vs_non_daps_AH46921 <- FT3(daps, non_daps, AH46921, "E068")
or_bp_daps_vs_non_daps_AH46922 <- FT3(daps, non_daps, AH46922, "E069")
or_bp_daps_vs_non_daps_AH46924 <- FT3(daps, non_daps, AH46924, "E071")
or_bp_daps_vs_non_daps_AH46926 <- FT3(daps, non_daps, AH46926, "E073")

or_bp_daps_vs_non_daps_chromHMM <- bind_rows(or_bp_daps_vs_non_daps_AH46921,
                                             or_bp_daps_vs_non_daps_AH46922,
                                             or_bp_daps_vs_non_daps_AH46924,
                                             or_bp_daps_vs_non_daps_AH46926) %>%
  mutate(feature = factor(feature,
                          levels = c("Active TSS", "Flanking Active TSS",
                                     "Transcr. at gene 5' and 3'",
                                     "Strong transcription",
                                     "Weak transcription", "Genic enhancers",
                                     "Enhancers", "ZNF genes & repeats",
                                     "Heterochromatin", "Bivalent/Poised TSS",
                                     "Flanking Bivalent TSS/Enh",
                                     "Bivalent Enhancer", "Repressed PolyComb",
                                     "Weak Repressed PolyComb",
                                     "Quiescent/Low"),
                          ordered = TRUE))

g <- or_bp_daps_vs_non_daps_chromHMM %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Base-level: DAPs vs. non-DAPs") +
  ylab("log2(OR) with 95% CI") +
  scale_colour_manual(values = cols)
ggsave("DAPs.chromHMM-enrichment.pdf", plot = g, width = 5, height = 5)


#-------------------------------------------------------------------------------
# DAPs vs. non-DAPs (directional): chromHMM 15-state model for 4 brain regions
#

or_bp_daps_vs_non_daps_AH46921_hypo <- cbind(
  FT3(daps_hypo, non_daps, AH46921, "E068"),
  data_frame(direction = "hypo"))
or_bp_daps_vs_non_daps_AH46921_hyper <- cbind(
  FT3(daps_hyper, non_daps, AH46921, "E068"),
  data_frame(direction = "hyper"))
or_bp_daps_vs_non_daps_AH46922_hypo <- cbind(
  FT3(daps_hypo, non_daps, AH46922, "E069"),
  data_frame(direction = "hypo"))
or_bp_daps_vs_non_daps_AH46922_hyper <- cbind(
  FT3(daps_hyper, non_daps, AH46922, "E069"),
  data_frame(direction = "hyper"))
or_bp_daps_vs_non_daps_AH46924_hypo <- cbind(
  FT3(daps_hypo, non_daps, AH46924, "E071"),
  data_frame(direction = "hypo"))
or_bp_daps_vs_non_daps_AH46924_hyper <- cbind(
  FT3(daps_hyper, non_daps, AH46924, "E071"),
  data_frame(direction = "hyper"))
or_bp_daps_vs_non_daps_AH46926_hypo <- cbind(
  FT3(daps_hypo, non_daps, AH46926, "E073"),
  data_frame(direction = "hypo"))
or_bp_daps_vs_non_daps_AH46926_hyper <- cbind(
  FT3(daps_hyper, non_daps, AH46924, "E073"),
  data_frame(direction = "hyper"))

or_bp_directional_daps_vs_non_daps_chromHMM <-
  bind_rows(or_bp_daps_vs_non_daps_AH46921_hypo,
            or_bp_daps_vs_non_daps_AH46921_hyper,
            or_bp_daps_vs_non_daps_AH46922_hypo,
            or_bp_daps_vs_non_daps_AH46922_hyper,
            or_bp_daps_vs_non_daps_AH46924_hypo,
            or_bp_daps_vs_non_daps_AH46924_hyper,
            or_bp_daps_vs_non_daps_AH46926_hypo,
            or_bp_daps_vs_non_daps_AH46926_hyper) %>%
  mutate(feature = factor(feature,
                          levels = c("Active TSS", "Flanking Active TSS",
                                     "Transcr. at gene 5' and 3'",
                                     "Strong transcription",
                                     "Weak transcription", "Genic enhancers",
                                     "Enhancers", "ZNF genes & repeats",
                                     "Heterochromatin", "Bivalent/Poised TSS",
                                     "Flanking Bivalent TSS/Enh",
                                     "Bivalent Enhancer", "Repressed PolyComb",
                                     "Weak Repressed PolyComb",
                                     "Quiescent/Low"),
                          ordered = TRUE))

g <- or_bp_directional_daps_vs_non_daps_chromHMM %>%
  ggplot(aes(x = feature, y = log2(estimate), col = db)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2(lower), ymax = log2(upper))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Base-level: DAPs vs. non-DAPs") +
  ylab("log2(OR) with 95% CI") +
  scale_colour_manual(values = cols) +
  facet_grid(~direction)

ggsave("DAPs-directional.chromHMM-enrichment.pdf", plot = g, width = 5,
       height = 5)
