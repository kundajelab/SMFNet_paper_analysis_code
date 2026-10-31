#Run from the repository root: Rscript figure_code/Fig2_E4_PRISM_cluster_benchmarking/plot_PRISM_validation.R
library(dplyr)
library(tidyr)
library(ggplot2)
library(egg)
library(gridExtra)

PLOT_DIR <- "output/Fig2_E4_PRISM_cluster_benchmarking/PRISM_validation"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
analysis_dir <- "input_files/Fig2_E4_PRISM_cluster_benchmarking"

load_real_agreement <- function(filename) {
  read.delim(file.path(analysis_dir, filename)) %>%
    filter(target_kind == "real", kernel == "fraction_agreement")
}

real_effects <- load_real_agreement("prism_effects.tsv")
real_bins <- load_real_agreement("prism_bins.tsv") %>%
  filter(selection_quantile %in% 0:2) %>%
  arrange(selection_quantile)
real_summary <- load_real_agreement("prism_summary.tsv") %>%
  dplyr::slice(1)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6),
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

#Panel A
cor_label <- sprintf(
  "rho=%.2f\n95%% CI [%.2f, %.2f]",
  real_summary$spearman,
  real_summary$spearman_ci_low,
  real_summary$spearman_ci_high
)

prism_scatter_val <- real_effects %>%
  ggplot(aes(x = prism_delta, y = 100 * delta_sim)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 0.3, color = "#454545", alpha = 0.3) +
  annotate("text", x = -Inf, y = Inf,
    label = cor_label,
    hjust = -0.1, vjust = 1.2,
    size = 5/2.8) +
  labs( x = "PRISM mutation score", y = "Change in methylation agreement (pp)")

pdf(paste0(PLOT_DIR, "/prism_validation_scatter.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(prism_scatter_val, height = unit(3, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
dev.off()

#Panel B
real_bins <- real_bins %>%
  mutate(
    group = factor(
      selection_quantile,
      levels = 0:2,
      labels = c("Most\nnegative", "Closest\nto zero", "Most\npositive")
    ),
    mean_pp = 100 * mean_delta_sim,
    lower_pp = 100 * delta_sim_ci_low,
    upper_pp = 100 * delta_sim_ci_high
  )

prism_sum_val <- real_bins %>%
  ggplot(aes(x = group, y = mean_pp, fill = group)) +
  geom_errorbar(aes(ymin = lower_pp, ymax = upper_pp), width = 0.3) +
  geom_point(shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#676767", "#BFBFBF", "#EFEFEF")) +
  labs(y = "Mean change in agreement (pp)") +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  coord_cartesian(ylim = c(-3.5, 3.5))

pdf(paste0(PLOT_DIR, "/prism_validation_sum.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(prism_sum_val, height = unit(2, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
dev.off()
