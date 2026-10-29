#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_SMFNet_pos_comp.R
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
  plot.title = element_text(size = 6),
  title = element_text(size = 6),
  axis.title = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank(),
  plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
))

#Expand the comma-separated site positions and average signals.
read_averages <- function(filename, header, source) {
  read.delim(file.path(input_dir, filename), header = header,
             col.names = c("region", "smf_pos", "avg_signal")) %>%
    separate_longer_delim(c(smf_pos, avg_signal), delim = ",") %>%
    mutate(type = source)
}
observed <- read_averages("avg_TestRegions_Experiment.tsv", FALSE, "act")
predicted <- read_averages("avg_TestRegions_SMFNet.tsv", TRUE, "pred")

#Match by region and site; retain the original scatter plot's complete pairs without coverage filtering.
averages <- bind_rows(observed, predicted) %>%
  pivot_wider(names_from = type, values_from = avg_signal) %>%
  na.omit() %>%
  mutate(act = as.numeric(act), pred = as.numeric(pred))
pearson <- cor(averages$act, averages$pred)

average_plot <- ggplot(averages, aes(x = act, y = pred)) +
  geom_hex(bins = 500) +
  scale_fill_viridis_c(option = "A", trans = "log10") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  annotate("text", x = 0.15, y = 0.95, size = 5 / 2.8,
           label = paste0("r=", round(pearson, 2))) +
  labs(x = "Observed average", y = "Predicted average")

pdf(file.path(output_dir, "SMFNet_vs_Experiment_smfpos_scatter.pdf"),
    width = 20 / 2.54, height = 20 / 2.54)
panel <- set_panel_size(average_plot, height = grid::unit(2, "cm"),
                        width = grid::unit(2, "cm"))
grid.arrange(panel)
dev.off()
