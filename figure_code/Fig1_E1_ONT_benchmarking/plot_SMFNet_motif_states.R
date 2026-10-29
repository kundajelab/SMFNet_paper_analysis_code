#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_SMFNet_motif_states.R
library(ggplot2)
library(viridis)
library(egg)
library(gridExtra)

#Use --correlations-only to plot state correlations without the motif-average table.
correlations_only <- "--correlations-only" %in% commandArgs(trailingOnly = TRUE)

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

save_panel <- function(plot, filename, height, width) {
  pdf(file.path(output_dir, filename), width = 20 / 2.54, height = 20 / 2.54)
  on.exit(dev.off())
  panel <- set_panel_size(plot, height = grid::unit(height, "cm"),
                          width = grid::unit(width, "cm"))
  grid.arrange(panel)
}

#Read the tables produced by count_SMFNet_motif_states.py and compute_SMFNet_motif_avg.py.
state_counts <- read.delim(file.path(output_dir, "SMFNet_motif_state_counts.tsv"))
state_counts <- state_counts[!is.na(state_counts$act_n_valid) & state_counts$act_n_valid > 100, ]

message("State-correlation motif instances after act_n_valid > 100: ", nrow(state_counts))

state_labels <- c(accessible = "ac", bound = "bound", nucleosomal = "nuc")
for (state in names(state_labels)) {
  observed <- paste0("act_", state, "_perc")
  predicted <- paste0("pred_", state, "_perc")
  pearson <- cor(state_counts[[observed]], state_counts[[predicted]])
  plot <- ggplot(state_counts, aes(x = .data[[observed]], y = .data[[predicted]])) +
    geom_hex(bins = 200) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_viridis(trans = "log10", option = "A") +
    annotate("text", label = paste0("r=", round(pearson, 2)),
             size = 5 / 2.8, x = 0.15, y = 0.9)
  save_panel(plot, paste0("cor_smfnet_", state_labels[[state]], ".pdf"), 2, 2)
}

if (!correlations_only) {
  composites <- read.delim(file.path(output_dir, "SMFNet_motif_composites_mean.tsv"))
  composite_plot <- ggplot(composites, aes(x = rel_pos, y = mean_value, color = dataset)) +
    facet_wrap(~motif) +
    geom_line() +
    scale_y_continuous(limits = c(0.2, 0.8)) +
    scale_color_manual(values = c("#BCBCBC", "#636cb2"))
  save_panel(composite_plot, "motif_composites.pdf", 1, 1)
}
