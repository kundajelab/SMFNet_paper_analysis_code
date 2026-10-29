#Run from SMFNet_paper_analysis_code github repo:
#Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_SMFNet_metrics.R
library(ggplot2)
library(dplyr)
library(tidyr)
library(egg)
library(gridExtra)

input_dir <- file.path("input_files", "Fig1_E1_ONT_benchmarking")
output_dir <- file.path("output", "Fig1_E1_ONT_benchmarking")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(),
  axis.text = element_text(size = 6),
  title = element_text(size = 6),
  axis.title = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank()
))

metrics <- read.delim(file.path(input_dir, "SMFNet_WT_SMF-ONT_metrics.tsv")) %>%
  filter(n_observed_reads >= 100)
message("Regions with at least 100 observed reads: ", nrow(metrics))

colors <- c(local_regression = "#BCBCBC", shuffle = "#BCBCBC",
            pseudorep = "#DEDEDE", smfnet_full = "#636bb2",
            smfnet_matched = "#9FA4D3")
labels <- c(local_regression = "Local\nregression", shuffle = "Shuffle",
            pseudorep = "Pseudo-\nreplicates", smfnet_full = "SMFNet\nfull",
            smfnet_matched = "SMFNet\nmatched")

#General purpose plotting function
save_violin <- function(data, comparisons, filename, ylabel, limits, width) {
  data <- data %>% filter(is.finite(value)) %>%
    mutate(comparison = factor(comparison, levels = comparisons))
  print(data %>% group_by(comparison) %>%
          summarise(n_regions = n(), median = median(value), .groups = "drop"))
  plot <- ggplot(data, aes(x = comparison, y = value, fill = comparison)) +
    geom_violin(color = NA) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = labels) +
    coord_cartesian(ylim = limits) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_summary(fun = median, geom = "crossbar", width = 0.3, linewidth = 0.25) +
    labs(x = NULL, y = ylabel) +
    theme(legend.position = "none")
  pdf(file.path(output_dir, filename), width = 20 / 2.54, height = 20 / 2.54)
  on.exit(dev.off())
  panel <- set_panel_size(plot, height = grid::unit(2.5, "cm"),
                          width = grid::unit(width, "cm"))
  grid.arrange(panel)
}

#Average profiles: retain the original local/flank-regression baseline column.
average_metrics <- metrics %>%
  transmute(region, local_regression = pearson_flank_baseline,
            pseudorep = pearson_matched_pseudo_replicate,
            smfnet_full = pearson, smfnet_matched = pearson_matched) %>%
  pivot_longer(-region, names_to = "comparison", values_to = "value")
save_violin(average_metrics,
            c("local_regression", "pseudorep", "smfnet_full", "smfnet_matched"),
            "TestRegions_average_profile_pearson.pdf", "Average-profile Pearson r",
            c(-0.4, 1), 2.5)

#Lengths of occupied runus (state 1): retain the original read-axis shuffle control.
length_metrics <- metrics %>%
  transmute(region, shuffle = run_length_ones_wasserstein_read_axis_shuffle,
            pseudorep = run_length_ones_wasserstein_pseudo_replicate,
            smfnet_matched = run_length_ones_wasserstein) %>%
  pivot_longer(-region, names_to = "comparison", values_to = "value")
save_violin(length_metrics, c("shuffle", "pseudorep", "smfnet_matched"),
            "SMFNet_mESC_ONT_TestRegions_occupied_length_wasserstein.pdf", "Occupied length Wasserstein distance",
            c(0, 55), 2)

#Pearson correlation between site-correlation matrices.
correlation_metrics <- metrics %>%
  transmute(region, pseudorep = covariation_corr_of_corr_pearson_pseudo_replicate,
            smfnet_matched = covariation_corr_of_corr_pearson,
            smfnet_full = covariation_corr_of_corr_pearson_all_reads) %>%
  pivot_longer(-region, names_to = "comparison", values_to = "value")
save_violin(correlation_metrics, c("pseudorep", "smfnet_matched", "smfnet_full"),
            "SMFNet_mESC_ONT_TestRegions_correlation_of_correlations.pdf", "Correlation of correlations (r)",
            c(0, 1), 2)
