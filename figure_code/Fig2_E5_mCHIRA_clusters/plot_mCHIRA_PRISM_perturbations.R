#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig2_E5_mCHIRA_clusters/plot_mCHIRA_PRISM_perturbations.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(egg)
library(gridExtra)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(),
                  plot.title = element_text(size = 6), axis.text = element_text(size = 6),
                  title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), plot.margin = grid::unit(c(0, 0, 0, 0), "cm")))


OUT_DIR <- "output/Fig2_E5_mCHIRA_clusters/"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#Compare motif PRISM scores with mutant-minus-wildtype cluster fractions.
new_cluster <- c("2"="A","1"="B","3"="C","0"="D","4"="E","5"="F","6"="G","7"="H")
prism_scores <- read.delim("output/Fig2_E5_mCHIRA_clusters/CRE_6005_cluster_PRISM_summary.tsv")
cluster_fractions <- read.delim("output/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_cluster_fractions.tsv")
cluster_fractions$cluster <- unname(new_cluster[as.character(cluster_fractions$cluster)])

wt_fractions <- cluster_fractions %>%
  dplyr::filter(group == "CRE_6005") %>%
  dplyr::transmute(sample, cluster, wt_fraction = fraction)

cluster_delta <- cluster_fractions %>%
  dplyr::left_join(wt_fractions) %>%
  dplyr::mutate(delta = fraction - wt_fraction)

nrf1_pos <- 103:112
myc_pos <- 178:185

nrf1_score <- prism_scores %>%
  dplyr::filter(position %in% nrf1_pos) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarize(prism = mean(signal)) %>%
  dplyr::mutate(group = "CRE_6005_mutant_1")

myc_score <- prism_scores %>%
  dplyr::filter(position %in% myc_pos) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarize(prism = mean(signal)) %>%
  dplyr::mutate(group = "CRE_6005_mutant_4")

scores <- rbind(nrf1_score, myc_score)

combi_df <- cluster_delta %>%
  dplyr::filter(group %in% c("CRE_6005_mutant_1", "CRE_6005_mutant_4")) %>%
  dplyr::left_join(scores) %>%
  dplyr::mutate(mutant = ifelse(group == "CRE_6005_mutant_1", "dNRF1", "dMAX::MYC"))

cor_df <- combi_df %>%
  dplyr::group_by(sample) %>%
  dplyr::summarize(cor_prism = cor(delta, prism, method = "pearson"))

prism_fracs <- combi_df %>%
  ggplot() +
  facet_wrap(~sample) +
  geom_point(aes(x = prism, y = delta, color = mutant), size = 0.8) +
  geom_text(data = cor_df, y = 0.1, x = 1.25, aes(label = paste0("r=", round(cor_prism, 2))), size = 5/2.8) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_x_continuous(limits = c(-0.4, 1.5)) +
  scale_y_continuous(limits = c(-0.4, 0.2)) +
  scale_color_manual(values = c("#abb6dd", "#056955"))

pdf(file.path(OUT_DIR,paste0("prism_fracs_cor_scatter.pdf")),
    width=20/2.54,height=20/2.54)

grid.arrange(set_panel_size(prism_fracs,height=grid::unit(2,"cm"),width=grid::unit(2,"cm")))
dev.off()
