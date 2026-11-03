#Run from the repository root: Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_TFBS_coaccessibility_scatter.R
#Run prepare_TFBS_pairs.py and compute_TFBS_coaccessibility.py first.
#The scatter uses background-corrected pairwise scores; Pearson r uses all complete pairs.
library(ggplot2)
library(egg)
library(gridExtra)

input_dir <- file.path("input_files", "Fig1_E3_Codependency_benchmarking")
output_dir <- file.path("output", "Fig1_E3_Codependency_benchmarking", "TFBS")
data <- read.delim(file.path(output_dir, "TFBS_coaccessibility.tsv"))
correction <- read.delim(file.path(input_dir, "coaccessibility_background_bias.txt"))
bias <- correction$bias[match(c("mi_act", "mi_pred"), correction$type)]
stopifnot(length(bias) == 2, all(is.finite(bias)))
data$mi_act <- data$mi_act - bias[1]
data$mi_pred <- data$mi_pred - bias[2]
correlation <- cor(data$mi_act, data$mi_pred, method = "pearson", use = "complete.obs")
message("Pearson r: ", correlation)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), plot.title = element_text(size = 6), axis.text = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6), strip.text = element_text(size = 6),
  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
))

#Keep the original 100-bin hexagonal density and displayed score range.
plot <- ggplot(data, aes(mi_act, mi_pred)) +
  geom_hex(bins = 100) +
  scale_x_continuous(limits = c(-0.4, 1), name = "Experiment co-accessibility") +
  scale_y_continuous(limits = c(-0.4, 1), name = "SMFNet co-accessibility") +
  scale_fill_viridis_c(option = "A", trans = "log10") + coord_equal() +
  annotate("text", label = paste0("r=", round(correlation, 2)), x = -0.1, y = 0.8, size = 6 / 2.8) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed")
path <- file.path(output_dir, "TFBS_coaccessibility_scatter.pdf")
pdf(path, width = 20 / 2.54, height = 20 / 2.54)
grid.arrange(set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm")))
dev.off()
message("Wrote: ", path)
