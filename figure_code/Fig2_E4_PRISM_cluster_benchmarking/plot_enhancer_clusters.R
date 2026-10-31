#Run from the repository root: Rscript figure_code/Fig2_E4_PRISM_cluster_benchmarking/plot_enhancer_clusters.R
library(dplyr)
library(ggplot2)
library(egg)
library(gridExtra)

PLOT_DIR <- "output/Fig2_E4_PRISM_cluster_benchmarking/"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(),
                  plot.title = element_text(size = 6), axis.text = element_text(size = 6),
                  title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), plot.margin = grid::unit(c(0, 0, 0, 0), "cm")))

sum_df <- read.delim("output/Fig2_E4_PRISM_cluster_benchmarking/enhancer_cluster_fractions.tsv")

cor_clusters <- cor(sum_df$pred_frac, sum_df$act_frac)

#Compare cluster fractions within each region.
region_metrics <- sum_df %>%
  group_by(region) %>%
  summarize(
    n_clusters = n(),
    r = cor(pred_frac, act_frac, use = "complete.obs"),
    .groups = "drop"
  )

median_r <- median(region_metrics$r)

cluster_pearson <- region_metrics %>%
  ggplot(aes(x = r)) +
  geom_histogram(bins = 60, fill = "#dedede") +
  labs(x = "Per-region Pearson", y = "Number of regions") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = median_r), linetype = "dashed", color = "#b2092d") +
  annotate(geom = "text", label = paste0(round(median_r, 2)), x = median_r - 0.2, y = 80, color = "#b2092d", size = 5/2.8)

pdf(paste0(PLOT_DIR, paste0("cluster_region_pearson.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(cluster_pearson, height = grid::unit(1.5, "cm"), width = grid::unit(2, "cm"))
grid.arrange(fix)
dev.off()

#Plot predicted minus experimental fractions against their mean.
sum_df_eval <- sum_df %>%
  mutate(
    mean_frac = (act_frac + pred_frac) / 2,
    diff_frac = pred_frac - act_frac
  )

diff_plot <- sum_df_eval %>%
  ggplot(aes(x = mean_frac, y = diff_frac)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hex(bins = 300) +
  scale_fill_viridis_c(trans = "log10", option = "A") +
  scale_x_continuous(name = "Mean fraction") +
  scale_y_continuous(name = "dSMFNet-Experiment (%)") +
  coord_cartesian(xlim = c(0, 0.6)) +
  annotate(geom = "text", label = paste0("r=", round(cor_clusters, 2)), size = 5/2.8, x = 0.1, y = 0.45)

pdf(paste0(PLOT_DIR, paste0("cluster_diff.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(diff_plot, height = grid::unit(3, "cm"), width = grid::unit(4.5, "cm"))
grid.arrange(fix)
dev.off()

