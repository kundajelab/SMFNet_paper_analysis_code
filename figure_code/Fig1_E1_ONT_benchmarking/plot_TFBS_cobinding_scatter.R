#Run from the repository root: Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_TFBS_cobinding_scatter.R
#Run prepare_TFBS_pairs.py and compute_TFBS_cobinding.py first.
#Odds ratios use bound state 2 and a 0.5 pseudocount.
library(ggplot2)
library(egg)
library(gridExtra)

output_dir <- file.path("output", "Fig1_E3_Codependency_benchmarking", "TFBS", "cobinding")
data <- read.delim(file.path(output_dir, "TFBS_cobinding.tsv"))
data <- data[as.logical(data$plot_eligible), ]
correlation <- cor(data$log_odds_act, data$log_odds_pred)
message("Pairs: ", nrow(data), "; Pearson r: ", correlation)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), axis.text = element_text(size = 6), plot.title = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6), strip.text = element_text(size = 6),
  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
))
plot <- ggplot(data, aes(log_odds_act, log_odds_pred)) +
  geom_hex(bins = 90) + scale_fill_viridis_c(option = "A", trans = "log10") +
  scale_x_continuous(limits = c(-3, 13), name = "Experiment log2 odds ratio") +
  scale_y_continuous(limits = c(-3, 13), name = "SMFNet log2 odds ratio") +
  annotate("text", x = 1, y = 10, size = 5 / 2.8, label = paste0("r=", round(correlation, 2))) +
  theme(legend.position = "right")
path <- file.path(output_dir, "TFBS_cobinding_log2_odds_scatter.pdf")
pdf(path, width = 20 / 2.54, height = 20 / 2.54)
grid.arrange(set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm")))
dev.off()
message("Wrote: ", path)
