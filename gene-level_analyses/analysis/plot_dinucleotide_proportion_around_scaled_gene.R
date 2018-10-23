# Plot dinucleotide proportions around genes using equal number of bins along
# gene body +/- 2 gene body equivalents (GBE)
# Peter Hickey
# 2018-03-07

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(cowplot)
library(broom)
library(tidyr)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "gene-level_analyses_data.rda"))
x <- readRDS(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "dinucleotide_count_around_around_scaled_gene.rds"))
y <- readRDS(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "dinucleotide_count_around_around_scaled_gene.opposite_strand.rds"))

### ----------------------------------------------------------------------------
### Data wrangling
###

sx <- x %>%
  inner_join(select(genes_df, gene, CGI_promoter, gene_type)) %>%
  group_by(dinuc, bin, CGI_promoter, gene_type) %>%
  summarise(proportion = weighted.mean(n_dinucleotide / width, width))

sx_logFC <- inner_join(
  x,
  select(genes_df, gene, CGI_promoter, logFC, gene_type) %>%
    mutate(qLogFC = cut(logFC,
                        quantile(logFC, 0:4 / 4, na.rm = TRUE),
                        ordered_result = TRUE))) %>%
  group_by(dinuc, bin, CGI_promoter, qLogFC, gene_type) %>%
  summarise(proportion = weighted.mean(n_dinucleotide / width, width))

sx_rpkm_NAcc_pos <- inner_join(
  x,
  select(genes_df, gene, CGI_promoter, rpkm_NAcc_pos, gene_type) %>%
    mutate(qrpkm_NAcc_pos = cut(rpkm_NAcc_pos,
                        quantile(rpkm_NAcc_pos, 0:4 / 4, na.rm = TRUE),
                        ordered_result = TRUE))) %>%
  group_by(dinuc, bin, CGI_promoter, qrpkm_NAcc_pos, gene_type) %>%
  summarise(proportion = weighted.mean(n_dinucleotide / width, width))

sx_rpkm_BA9_pos <- inner_join(
  x,
  select(genes_df, gene, CGI_promoter, rpkm_BA9_pos, gene_type) %>%
    mutate(qrpkm_BA9_pos = cut(rpkm_BA9_pos,
                                quantile(rpkm_BA9_pos, 0:4 / 4, na.rm = TRUE),
                                ordered_result = TRUE))) %>%
  group_by(dinuc, bin, CGI_promoter, qrpkm_BA9_pos, gene_type) %>%
  summarise(proportion = weighted.mean(n_dinucleotide / width, width))

sy <- y %>%
  inner_join(select(genes_df, gene, CGI_promoter, gene_type)) %>%
  group_by(dinuc, bin, CGI_promoter, gene_type) %>%
  summarise(proportion = weighted.mean(n_dinucleotide / width, width))

### ----------------------------------------------------------------------------
### Plots
###

g <- ggplot(sx, aes(x = bin, y = proportion, col = CGI_promoter)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 200, lty = 2) +
  geom_vline(xintercept = 300, lty = 2) +
  geom_line() +
  facet_wrap(gene_type ~ dinuc, ncol = 8) +
  scale_colour_brewer(palette = "Dark2") +
  ylim(0, NA) +
  ggtitle("Same strand")
save_plot("../figures/dinucleotide_proportion_around_scaled_gene.pdf",
          g,
          base_height = 10,
          base_width = 10)

g <- ggplot(filter(sx_logFC, grepl("^C", dinuc), !is.na(qLogFC)),
            aes(x = bin, y = proportion, col = qLogFC)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 200, lty = 2) +
  geom_vline(xintercept = 300, lty = 2) +
  geom_line() +
  facet_wrap(gene_type + CGI_promoter ~ dinuc, ncol = 4,
             labeller = labeller(CGI_promoter = label_both)) +
  scale_colour_brewer(palette = "Dark2") +
  ylim(0, NA) +
  ggtitle("Same strand")
save_plot("../figures/dinucleotide_proportion_around_scaled_gene.logFC.pdf",
          g,
          base_height = 10,
          base_width = 10)

g <- ggplot(filter(sx_rpkm_NAcc_pos, grepl("^C", dinuc), !is.na(qrpkm_NAcc_pos)),
            aes(x = bin, y = proportion, col = qrpkm_NAcc_pos)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 200, lty = 2) +
  geom_vline(xintercept = 300, lty = 2) +
  geom_line() +
  facet_wrap(gene_type + CGI_promoter ~ dinuc, ncol = 4,
             labeller = labeller(CGI_promoter = label_both)) +
  scale_colour_brewer(palette = "Dark2") +
  ylim(0, NA) +
  ggtitle("Same strand")
save_plot("../figures/dinucleotide_proportion_around_scaled_gene.NAcc_pos_rpkm.pdf",
          g,
          base_height = 10,
          base_width = 10)

g <- ggplot(filter(sx_rpkm_BA9_pos, grepl("^C", dinuc), !is.na(qrpkm_BA9_pos)),
            aes(x = bin, y = proportion, col = qrpkm_BA9_pos)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 200, lty = 2) +
  geom_vline(xintercept = 300, lty = 2) +
  geom_line() +
  facet_wrap(gene_type + CGI_promoter ~ dinuc, ncol = 4,
             labeller = labeller(CGI_promoter = label_both)) +
  scale_colour_brewer(palette = "Dark2") +
  ylim(0, NA) +
  ggtitle("Same strand")
save_plot("../figures/dinucleotide_proportion_around_scaled_gene.BA9_pos_rpkm.pdf",
          g,
          base_height = 10,
          base_width = 10)

g <- ggplot(sy, aes(x = bin, y = proportion, col = CGI_promoter)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 200, lty = 2) +
  geom_vline(xintercept = 300, lty = 2) +
  geom_line() +
  facet_wrap(gene_type ~ dinuc, ncol = 8) +
  scale_colour_brewer(palette = "Dark2") +
  ylim(0, NA) +
  ggtitle("Opposite strand")
save_plot(
  "../figures/dinucleotide_proportion_around_scaled_gene.opposite_strand.pdf",
  g,
  base_height = 10,
  base_width = 10)
