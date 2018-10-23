# Plot LDSC results for 'adjusting for baseline' analyses
# Peter Hickey
# 2018-08-07

### ----------------------------------------------------------------------------
### Packages
###

library(GenomicRanges)
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
# NOTE: Only using a subset of features (e.g., not using bigDARs)
interesting_features <- c("POS_CG-DMRs", "POSvsNEG_CG-DMRs",
                          "CH-DMRs",
                          "POS_DARs", "POSvsNEG_DARs",
                          "H3K27ac", "chromHMM_union", "CNS")
brain_categories <- brain_categories[interesting_features]
brain_categories_df <- read_csv("../tables/brain_categories_df.csv") %>%
  mutate(Source = factor(Source, levels = c("NeuN+", "NeuN-", "NeuN+ & NeuN-",
                                            "Other", "Baseline")))

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
    x <- read_tsv(paste0("../output/ldsc/",
                         bc,
                         ".Height.Phase1.results")) %>%
      filter(Category == "L2_0" | Category == "CNS_0")
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
  fls <- list.files("../output/ldsc",
             pattern = glob2rx(
               paste0(bc, ".*Phase1.results")),
             full.names = TRUE)
  grep("adjusting", fls, invert = TRUE, value = TRUE)
}))

# Read in files, tidy up, and rbind
x <- bind_rows(lapply(fls, function(fl) {
  suppressMessages(read_tsv(fl)) %>%
    filter(Category == "L2_0" | Category == "CNS_0") %>%
    mutate(Category = sapply(strsplit(basename(fl), "\\."), "[[", 1),
           Trait = sapply(strsplit(sub(".Phase1.results", "", basename(fl)),
                                   "\\."),
                          "[[", 2),
           lower = Enrichment - 1.96 * Enrichment_std_error,
           upper = Enrichment + 1.96 * Enrichment_std_error,
           file = fl)
}))

# Join munged LDSC output with categories_df
x <- x %>%
  mutate(Category = factor(Category, levels(categories_df$Category),
                           ordered = TRUE)) %>%
  inner_join(categories_df, by = c("Category" = "Category")) %>%
  inner_join(traits_df, by = c("Trait" = "Trait"))

stopifnot(length(fls) == nrow(x))

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
  group_by(Category, file) %>%
  filter(grepl(Category, file)) %>%
  ungroup() %>%
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
###

strata_df <- x %>%
  group_by(Trait) %>%
  summarise(strata = ifelse(any(Coefficient_holm < 0.05),
                            "Brain-linked (sig)",
                            "Brain-linked (non-sig)"),
            TraitType2 = unique(TraitType2),
            strata = ifelse(TraitType2 == "Additional_phenotype",
                            "Non-brain-linked",
                            strata)) %>%
  dplyr::select(Trait, strata)
saveRDS(strata_df, "../objects/trait_strata_df.rds")
x_stratified <- inner_join(x, strata_df)

### ----------------------------------------------------------------------------
### Tables of results
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
  write_csv("../tables/LDSC_results.baseline_adjustments.csv")

x_stratified %>%
  dplyr::select(`Pretty Trait`, TraitType, N, N_cases, N_controls) %>%
  distinct() %>%
  rename(Trait = `Pretty Trait`,
         `Type` = `TraitType`) %>%
  write_csv("../tables/annotated_traits_df.csv")

### ----------------------------------------------------------------------------
### Plots
###

# ------------------------------------------------------------------------------
# Coefficient Z-score
#

# NOTE: ylim for coefficient Z-score plots chosen to fit full range of data
#       across baseline_adjustments, stringent_adjustments, and
#       adjusting_for_ndf_excluding_df
ylim_coefficient_score <- c(-4.5, 9)
g <- x_stratified %>%
  arrange(`Pretty Category`) %>%
  ggplot(data = .,
         aes(x = `Pretty Category`, y = `Coefficient_z-score`,
             col = `Pretty Category`, shape = sig_coef, size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ `Pretty Trait`, ncol = 5) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  scale_color_brewer(palette = "Dark2") +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  ylim(ylim_coefficient_score)
ggsave("../figures/Coefficient_Z-score.baseline_adjustments.pdf",
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
ggsave("../figures/Coefficient_Z-score.baseline_adjustments.stratified.pdf",
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
  arrange(`Pretty Category`) %>%
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
ggsave("../figures/Enrichment.baseline_adjustments.pdf",
       g,
       height = 6,
       width = 7)

# NOTE: ylim for this enrichment plot chosen from baseline_adjustments data
#       (i.e. no need to manually specify it)
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
  guides(col = FALSE, shape = FALSE, size = FALSE)
ggsave("../figures/Enrichment.baseline_adjustments.sig_stratified.pdf",
       g,
       height = 4,
       width = 5)

### ----------------------------------------------------------------------------
### Create legend used in all plot_ldsc.* output
###

g <- x %>%
  arrange(`Pretty Category`) %>%
  ggplot(data = .,
         aes(x = `Pretty Category`, y = -log10(Coefficient_p),
             col = `Pretty Category`, shape = sig_coef, size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ Trait, ncol = 5) +
  # Holm's cutoff
  geom_hline(aes(yintercept = -log10(Coefficient_holm_cutoff))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3))
legend_plot <- ggdraw(plot_grid(NULL, get_legend(g)))
ggsave("../figures/Legend.pdf",
       legend_plot,
       height = 6,
       width = 6)
