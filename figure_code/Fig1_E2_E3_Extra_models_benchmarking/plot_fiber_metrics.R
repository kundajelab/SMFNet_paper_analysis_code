#Run from the repository root: Rscript figure_code/Fig1_E2_E3_Extra_models_benchmarking/plot_fiber_metrics.R
library(ggplot2)
library(dplyr)
library(tidyr)
library(egg)
library(gridExtra)

input_dir <- file.path("input_files", "Fig1_E2_E3_Extra_models_benchmarking")
output_dir <- file.path("output", "Fig1_E2_E3_Extra_models_benchmarking")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), axis.text = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6),
  strip.text = element_text(size = 6), strip.background = element_blank()
))

#Load the Fiber-seq model metrics independently of the Kreibich models.
fiber_metrics <- read.delim(file.path(input_dir, "GM12878_Fiber-seq_SMFNet_metrics.tsv")) %>%
  mutate(cell_type = "GM12878")

#Keep the original complete-case cohorts, including checks on unplotted controls.
#Those columns affect region selection but are not included in the plot data.
prepare_metrics <- function(data, metric) {
  if (metric == "average") {
    columns <- c(pearson_flanks = "pearson_flank_baseline",
                 pearson_pseudoreps = "pearson_matched_pseudo_replicate",
                 pearson = "pearson", pearson_matched = "pearson_matched")
    required <- c(unname(columns), "pearson_positional_baseline")
    comparisons <- names(columns)
  } else {
    columns <- c(covar_pseudoreps = "covariation_corr_of_corr_pearson_pseudo_replicate",
                 covar_all = "covariation_corr_of_corr_pearson_all_reads",
                 covar_matched = "covariation_corr_of_corr_pearson")
    required <- c(unname(columns), "covariation_corr_of_corr_pearson_positional_baseline")
    comparisons <- names(columns)
  }
  data %>%
    filter(n_observed_reads >= 100) %>%
    drop_na(all_of(c("region", required))) %>%
    select(region, cell_type, all_of(columns[comparisons])) %>%
    pivot_longer(-c(region, cell_type), names_to = "comparison", values_to = "value") %>%
    mutate(comparison = factor(comparison, levels = comparisons))
}

colors <- c(pearson_flanks = "#BCBCBC", pearson_pseudoreps = "#DEDEDE",
            pearson = "#636bb2", pearson_matched = "#9FA4D3",
            covar_pseudoreps = "#DEDEDE", covar_all = "#636bb2", covar_matched = "#9FA4D3")

#Use identical panel dimensions and median markers across model comparisons.
save_violin <- function(data, filename, ylabel, limits, width) {
  print(data %>% group_by(cell_type, comparison) %>%
          summarise(n_regions = n(), median = median(value), .groups = "drop"))
  plot <- ggplot(data, aes(x = comparison, y = value, fill = comparison)) +
    geom_violin(color = NA) +
    scale_fill_manual(values = colors) +
    coord_cartesian(ylim = limits) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_summary(fun = median, geom = "crossbar", width = 0.3, linewidth = 0.25) +
    labs(x = NULL, y = ylabel) +
    theme(legend.position = "none")
  pdf(file.path(output_dir, filename), width = 20 / 2.54, height = 20 / 2.54)
  on.exit(dev.off())
  panel <- set_panel_size(plot, height = grid::unit(2, "cm"), width = grid::unit(width, "cm"))
  grid.arrange(panel)
}

#Show both full and read-count-matched predictions for Fiber-seq.
save_violin(prepare_metrics(fiber_metrics, "average"),
            "TestRegions_avg_pearson_fiber.pdf", "Average-profile Pearson r",
            c(0, 1), 2)
save_violin(prepare_metrics(fiber_metrics, "covariation"),
            "TestRegions_avg_covar_fiber.pdf", "Correlation of correlations (r)",
            c(0, 1), 1.5)
