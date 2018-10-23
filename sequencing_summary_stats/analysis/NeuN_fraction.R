# Investigate how PMI and age of sample influences NeuN+ fraction in FANS
# Peter Hickey
# 2018-03-20

library(dplyr)
library(ggplot2)
library(lubridate)
library(lme4)
library(cowplot)

# ------------------------------------------------------------------------------
# Load data
#

# Read in data and convert some variables into a more friendly form
x <- read.csv("../extdata/NICHD_POS_proportions_from_notebook.csv",
              stringsAsFactors = FALSE) %>%
  mutate(Individual = factor(Individual),
         Tissue = gsub("Nacc", "NAcc", Tissue),
         Tissue = factor(Tissue),
         Date = ymd(Date),
         Age = sapply(strsplit(Age, ":"), function(x) {
           as.integer(x[[1]]) + as.integer(x[[2]]) / 365
         }),
         Sex = factor(Sex)) %>%
  as_data_frame() %>%
  arrange(Individual, Tissue) %>%
  group_by(Individual, Tissue) %>%
  mutate(Extraction = row_number()) %>%
  ungroup() %>%
  mutate(Extraction = factor(Extraction))

# ------------------------------------------------------------------------------
# Plots
#

# Plot the NeuN+ proportion as a function of (Individual, Age, PMI)
# stratified by Tissue

# Individual
g <- ggplot(x, aes(x = Tissue, y = POS / TOTAL, col = Individual)) +
  geom_jitter(width = 0.1) +
  # stat_summary(fun.data = "mean_cl_boot", colour = "black", geom = "crossbar") +
  ylab("NeuN+ Nuclei / Total Nuclei") +
  xlab("Brain Region")
save_plot("../figures/NeuN_fraction.pdf",
          g,
          base_height = 7)

# Age
g <- ggplot(x, aes(x = Age, y = POS / TOTAL)) +
  facet_wrap(~ Tissue) +
  geom_point(aes(col = Individual)) +
  geom_smooth(method = "lm") +
  # TODO: Why isn't SE showing up?
  # geom_smooth(
  #   method = "glm",
  #   aes(weight = TOTAL),
  #   method.args = list(family = "binomial"),
  #   se = TRUE) +
  ylab("NeuN Nuclei / Total Nuclei") +
  xlab("Age (years)")
save_plot("../figures/NeuN_fraction_by_age.pdf",
          g,
          base_height = 7)

# PMI
g <- ggplot(x, aes(x = PMI, y = POS / TOTAL)) +
  facet_wrap(~ Tissue) +
  geom_point(aes(col = Individual)) +
  geom_smooth(method = "lm") +
  # TODO: Why isn't SE showing up?
  # geom_smooth(
  #   method = "glm",
  #   aes(weight = TOTAL),
  #   method.args = list(family = "binomial"),
  #   se = TRUE) +
  ylab("NeuN Nuclei / Total Nuclei") +
  xlab("Post-mortem interval (hours)")
save_plot("../figures/NeuN_fraction_by_PMI.pdf",
          g,
          base_height = 7)

# ==============================================================================
# LMs/GLMs/GAMs
#

# ------------------------------------------------------------------------------
# Fitting LMs/GLMs/LMMs/GLMMs
#

# NOTE: All proportions are estimated with very high precision since the
#       min(TOTAL) is 1,344,887. Furthemore, the proportions are almost all away
#       from zero. Therefore, it is probably safe to use linear regression as an
#       approximation to the logistic regression (with the bonus of easier
#       interpretability).

# NOTE: lm() and glm() ignore the repeated measures whereas lmer() and glmer()
#       take this into account.
lm_fit <- lm(
  formula = POS / TOTAL ~ Tissue + PMI + Age,
  data = x)
lmer_fit <- lmer(
  formula = POS / TOTAL ~ Tissue + PMI + Age + (1 | Individual),
  data = x)
glm_fit <- glm(
  formula = cbind(POS, NEG) ~ Tissue + PMI + Age,
  data = x,
  family = binomial())
glmer_fit <- glmer(
  formula = cbind(POS, NEG) ~ Tissue + PMI + Age + (1 | Individual),
  data = x,
  family = binomial())

# ------------------------------------------------------------------------------
# Understanding GLM
#

# Intercept-only
intercept_only <- glm(cbind(POS, NEG) ~ 1, data = x, family = binomial())
# log(p / (1 - p)) = log(p) - log(1 - p) = logit(p) is the coefficient of this model
# Which means that the baseline odds are given by:
exp(coef(intercept_only))
# Which can be calculated from the raw data as:
p <- summarise(x, POS = sum(POS), NEG = sum(NEG)) %>%
  summarise(p = POS / (POS + NEG)) %>%
  pull(p)
p / (1 - p)

# Tissue-only
tissue_only <- glm(cbind(POS, NEG) ~ 0 + Tissue, data = x, family = binomial())
exp(coef(tissue_only))
p <- summarise(group_by(x, Tissue), POS = sum(POS), NEG = sum(NEG)) %>%
  mutate(p = POS / (POS + NEG)) %>%
  pull(p)
p / (1 - p)

# NOTE: These two models are identical!
#       `a` implictly aggregates data from the replicates by summing `POS`, `NEG
#       `b` explicitly aggregates data from the replicates by summing `POS`, `NEG
a <- glm(
  formula = cbind(POS, NEG) ~ Tissue * PMI,
  data = x,
  family = binomial())
b <- glm(
  formula = cbind(POS, NEG) ~ Tissue * PMI,
  data = x %>%
    group_by(Tissue, Individual, PMI) %>%
    summarise(POS = sum(POS),
              NEG = sum(NEG)),
  family = binomial())
