# Plot LDSC results for 'adjusting_for_ndf_excluding_df' analyses
# Peter Hickey
# 2018-08-07

### ----------------------------------------------------------------------------
### Packages
###

library(GenomicRanges)
library(rtracklayer)
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(gplots)
library(cowplot)

load("../../integrating_dmrs_dars_and_degs/objects/assays-and-features.rda")

###-----------------------------------------------------------------------------
### Category sizes
###

# ------------------------------------------------------------------------------
# Brain categories
#

brain_categories <- readRDS("../objects/categories.rds")
ndf_excluding_df <- readRDS("../objects/non-differential_features_excluding_differential_features.rds")
# NOTE: Only using a subset of features (e.g., not using bigDARs or
#       non-differential features)
interesting_features <- c("POS_CG-DMRs", "POSvsNEG_CG-DMRs",
                          "CH-DMRs",
                          "POS_DARs", "POSvsNEG_DARs")
brain_categories <- brain_categories[interesting_features]
brain_categories_df <- read_csv("../tables/brain_categories_df.csv") %>%
  mutate(Source = factor(Source, levels = c("NeuN+", "NeuN-", "NeuN+ & NeuN-",
                                            "Other", "Baseline"))) %>%
  filter(Category %in% names(brain_categories))

# ------------------------------------------------------------------------------
# Record proportion of SNPs and CpGs in each category
#

categories <- brain_categories

categories_df <- bind_rows(
  lapply(names(categories), function(x) {
    data_frame(Category = x,
               `Total width (bp)` = sum(width(categories[[x]])),
               `Mean width (bp)` = mean(width(categories[[x]])),
               `Median width (bp)` = median(width(categories[[x]])),
               n = length(categories[[x]]))
  })) %>%
  arrange(`Total width (bp)`) %>%
  inner_join(brain_categories_df)

snp_prop_table <- bind_rows(
  lapply(names(brain_categories), function(bc) {
    if (bc == "CNS") {
      x <- read_tsv("../output/ldsc/CNS.Height.Phase1.results",
                    n_max = 1)
    } else {
      x <- read_tsv(paste0("../output/ldsc/",
                           bc,
                           ".Height.Phase1.results"),
                    n_max = 1)
    }
    data_frame(Category = bc,
               `Prop._SNPs` = unlist(x[, "Prop._SNPs"]))
  }))
colnames(snp_prop_table)[2] <- "Proportion of SNPs"

cpg_prop_table <-
  data_frame(Category = names(categories),
             `Proportion of CpGs` = sapply(categories, function(x) {
               sum(overlapsAny(cpgs, x)) /
                 length(cpgs)}))

categories_df <- categories_df %>%
  inner_join(snp_prop_table, c("Category" = "Category")) %>%
  inner_join(cpg_prop_table, c("Category" = "Category")) %>%
  arrange(-Differential, `Total width (bp)`) %>%
  mutate(Category = factor(Category, Category, ordered = TRUE),
         `Pretty Category` = factor(`Pretty Category`, `Pretty Category`,
                                    ordered = TRUE)) %>%
  arrange(`Pretty Category`)

traits_df <- read_csv("../tables/traits_df.csv") %>%
  mutate(N = ifelse(is.na(N_cases), N, N_cases + N_controls),
         TraitType2 = ifelse(Trait == "BMI", "BMI", TraitType)) %>%
  mutate(TraitColour = case_when(
    .$TraitType2 == "Additional_phenotype" ~ brewer_pal("qual")(5)[1],
    .$TraitType2 == "Behavioural-cognitive" ~ brewer_pal("qual")(5)[2],
    .$TraitType2 == "Neurological" ~ brewer_pal("qual")(5)[3],
    .$TraitType2 == "Psychiatric" ~ brewer_pal("qual")(5)[4]),
    Trait2Colour = case_when(
      .$TraitType2 == "Additional_phenotype" ~ brewer_pal("qual")(5)[1],
      .$TraitType2 == "Behavioural-cognitive" ~ brewer_pal("qual")(5)[2],
      .$TraitType2 == "Neurological" ~ brewer_pal("qual")(5)[3],
      .$TraitType2 == "Psychiatric" ~ brewer_pal("qual")(5)[4],
      .$TraitType2 == "BMI" ~ brewer_pal("qual")(5)[5]))

categories_df <- categories_df %>%
  mutate(Group = case_when(.$Source == "NeuN+" & .$Differential &
                             !grepl("block", .$Category) ~ 1L,
                           .$Differential & !grepl("block", .$Category) ~ 2L,
                           .$Source == "NeuN+" & .$Differential ~ 3L,
                           .$Differential ~ 4L,
                           .$Source == "Baseline" ~ 6L,
                           TRUE ~ 5L),
         GroupColour = case_when(
           Group == 1L ~ "darkgreen",
           Group == 2L ~ "dodgerBlue",
           Group == 3L ~ "darkgreen",
           Group == 4L ~ "dodgerBlue",
           Group == 5L ~ "darkblue",
           Group == 6L ~ "darkgrey"))

### ----------------------------------------------------------------------------
### Load data and construct objects
###

fls <- unlist(lapply(names(brain_categories), function(bc) {
  list.files("../output/ldsc",
             pattern = glob2rx(
               paste0(bc, ".adjusting_for_ndf_excluding_df.*Phase1.results")),
             full.names = TRUE)
}))

# Read in files, tidy up, and rbind
x <- bind_rows(lapply(fls, function(fl) {
  suppressMessages(read_tsv(fl)) %>%
    filter(Category %in% c("L2_1", "L2_2", "L2_3", "L2_4")) %>%
    mutate(Category = case_when(
      .$Category == "L2_4" ~ sapply(strsplit(basename(fl), "\\."), "[[", 1),
      .$Category == "L2_1" ~ "CNS",
      .$Category == "L2_2" ~ "chromHMM_union",
      .$Category == "L2_3" ~ "H3K27ac",
      TRUE ~ gsub("_0", "", Category)),
      Trait = sapply(strsplit(sub(".Phase1.results", "", basename(fl)),
                              "\\."),
                     "[[", 3),
      lower = Enrichment - 1.96 * Enrichment_std_error,
      upper = Enrichment + 1.96 * Enrichment_std_error,
      file = fl)
}))

# Join munged LDSC output with categories_df
x <- x %>%
  mutate(Category = factor(Category,
                           levels(categories_df$Category),
                           ordered = TRUE)) %>%
  inner_join(categories_df, by = c("Category" = "Category")) %>%
  inner_join(traits_df, by = c("Trait" = "Trait"))

# NOTE: Anttila report these traits "had in sufficient evidence of additive
#       heritability for robust analysis" and excluded them from further
#       analysis
x <- x %>%
  filter(Trait != "Agreeableness",
         Trait != "Cardioembolic_stroke",
         Trait != "Large-vessel_disease",
         Trait != "Small-vessel_disease")

# Add adjusted P-values
x <- x %>%
  filter(Source != "Baseline") %>%
  mutate(Coefficient_p = pnorm(`Coefficient_z-score`, lower.tail = FALSE)) %>%
  group_by(Trait) %>%
  mutate(Coefficient_holm = p.adjust(Coefficient_p, method = "holm"),
         Enrichment_holm = p.adjust(Enrichment_p, method = "holm"),
         Coefficient_holm_cutoff =
           max(Coefficient_p[Coefficient_holm < 0.05], na.rm = TRUE),
         Enrichment_holm_cutoff =
           max(Enrichment_p[Enrichment_holm < 0.05], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sig_coef = Coefficient_holm < 0.05) %>%
  arrange(-Differential, `Total width (bp)`) %>%
  mutate(`Pretty Category` = factor(`Pretty Category`,
                                    unique(`Pretty Category`),
                                    ordered = TRUE),
         `Pretty Trait` = gsub("_", " ", Trait)) %>%
  arrange(`Pretty Category`)

### ----------------------------------------------------------------------------
### Stratify traits by 'Brain-linked (sig)', 'Brain-linked (non-sig)', or
### 'Non-brain-linked' (stratification defined in 'baseline' analysis)

strata_df <- readRDS("../objects/trait_strata_df.rds")
x_stratified <- inner_join(x, strata_df)

### ----------------------------------------------------------------------------
### Table of results
###

x_stratified %>%
  dplyr::select(-Category, -lower, -upper, -n, -`Prop._SNPs`,
                -file, -Source, -Differential, -Group, -GroupColour,
                -TraitType2, -TraitColour, -Trait2Colour, -Enrichment_p,
                -Enrichment_holm, -Enrichment_holm_cutoff,
                -`Mean width (bp)`, -`Median width (bp)`,
                -`Proportion of CpGs`,
                -Coefficient_holm_cutoff, -sig_coef, -Trait,
                -N, -N_cases, -N_controls, -TraitType) %>%
  dplyr::select(`Pretty Trait`, strata,
                `Pretty Category`, `Total width (bp)`, `Proportion of SNPs`,
                starts_with("Coefficient"), everything()) %>%
  rename(Feature = `Pretty Category`,
         Trait = `Pretty Trait`,
         `Proportion of h2` = `Prop._h2`,
         `Proportion of h2 standard error` = `Prop._h2_std_error`,
         Stratum = `strata`) %>%
  write_csv("../tables/LDSC_results.adjusting_for_ndf_excluding_df.csv")

### ----------------------------------------------------------------------------
### Plots
###

# NOTE: ylim for coefficient Z-score plots chosen to fit full range of data
#       across baseline_adjustments, stringent_adjustments, and
#       adjusting_for_ndf_excluding_df
ylim_coefficient_score <- c(-4.5, 9)

g <- x_stratified %>%
  filter(strata == "Brain-linked (sig)") %>%
  arrange(`Pretty Category`) %>%
  ggplot(data = .,
         aes(x = `Pretty Category`, y = `Coefficient_z-score`,
             col = `Pretty Category`, shape = sig_coef, size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ `Pretty Trait`, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_color_brewer(palette = "Dark2") +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  ylim(ylim_coefficient_score)
ggsave("../figures/Coefficient_Z-score.adjusting_for_ndf_excluding_df.pdf",
       g,
       height = 6,
       width = 7)

g <- x_stratified %>%
  ggplot(data = .,
         aes(x = `Pretty Category`, y = `Coefficient_z-score`,
             col = `Pretty Category`, shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(. ~ strata) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_color_brewer(palette = "Dark2") +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  ylim(ylim_coefficient_score)
ggsave("../figures/Coefficient_Z-score.adjusting_for_ndf_excluding_df.stratified.pdf",
       g,
       height = 4,
       width = 5)

# ------------------------------------------------------------------------------
# Enrichment
#

# NOTE: ylim for this enrichment plot chosen to fit full range of data across
#       baseline_adjustments, stringent_adjustments, and
#       adjusting_for_ndf_excluding_df
ylim_enrichment <- c(-7, 45)
g <- x_stratified %>%
  filter(strata == "Brain-linked (sig)") %>%
  ggplot(data = .,
         aes(x = `Pretty Category`, y = Enrichment, col = `Pretty Category`,
             shape = sig_coef)) +
  geom_point() +
  geom_pointrange(aes(ymin = Enrichment - 2 * Enrichment_std_error,
                      ymax = Enrichment + 2 * Enrichment_std_error)) +
  facet_wrap( ~ `Pretty Trait`, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, lty = 2) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  coord_cartesian(ylim = ylim_enrichment)
ggsave("../figures/Enrichment.adjusting_for_ndf_excluding_df.pdf",
       g,
       height = 6,
       width = 7)

# NOTE: ylim for this enrichment plot chosen to fit full range of data across
#       stringent_adjustments and adjusting_for_ndf_excluding_df
ylim_enrichment <- c(-7, 32)
g <- x_stratified %>%
  ggplot(data = .,
         aes(x = `Pretty Category`, y = Enrichment, col = `Pretty Category`,
             shape = sig_coef, size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ strata, labeller = labeller(sig = label_both)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_color_brewer(palette = "Dark2") +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  ylim(ylim_enrichment)
ggsave("../figures/Enrichment.adjusting_for_ndf_excluding_df.sig_stratified.pdf",
       g,
       height = 4,
       width = 5)
