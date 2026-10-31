#Run from the repository root: Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_TFBS_average_profiles.R
#Run prepare_TFBS_pairs.py and compute_TFBS_coaccessibility.py first.
#Uses the original raw mean SMF profiles.
library(ggplot2)
library(dplyr)
library(tidyr)
library(egg)
library(gridExtra)

input_dir <- file.path("input_files", "Fig1_E3_Codependency_benchmarking")
analysis_dir <- file.path("output", "Fig1_E3_Codependency_benchmarking", "TFBS")
output_dir <- file.path(analysis_dir, "average_profiles")
co_ac <- read.delim(file.path(analysis_dir, "TFBS_coaccessibility.tsv"))
pred_avg <- read.delim(file.path(input_dir, "avg_TestRegions_SMFNet.tsv"))
act_avg <- read.delim(file.path(input_dir, "avg_TestRegions_Experiment.tsv"), header = FALSE,
                      col.names = c("region", "smf_pos", "avg_signal"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), plot.title = element_text(size = 6), axis.text = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6), strip.text = element_text(size = 6),
  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
))

motifs_for_region <- function(region) {
  pairs <- co_ac[co_ac$region == region, ]
  bind_rows(transmute(pairs, motif = motif1, pos = pos1),
            transmute(pairs, motif = motif2, pos = pos2)) %>% distinct()
}

mean_profile <- function(data, region_id, label) {
  data[data$region == region_id, ] %>%
    separate_longer_delim(c(smf_pos, avg_signal), delim = ",") %>%
    mutate(smf_pos = as.integer(smf_pos), avg_signal = as.numeric(avg_signal), type = label)
}

save_average <- function(region, filename, offset = 0, width = 1024, motif_half_width = 5) {
  motifs <- motifs_for_region(region) %>% mutate(pos = pos - offset)
  data <- bind_rows(mean_profile(act_avg, region, "act"), mean_profile(pred_avg, region, "pred")) %>%
    mutate(smf_pos = smf_pos - offset)
  #Keep the original cropped example's inclusive 0..512 bp selection.
  if (width == 512) data <- data %>% filter(between(smf_pos, 0, width))
  plot <- ggplot(data) +
    geom_line(aes(smf_pos, avg_signal, color = type)) +
    geom_point(aes(smf_pos, avg_signal, color = type), size = 0.5) +
    geom_rect(data = motifs, aes(xmin = pos - motif_half_width, xmax = pos + motif_half_width,
      ymin = 0, ymax = 0.1), fill = "#BCBCBC") +
    geom_text(data = motifs, aes(pos, 0.15, label = motif), size = 4 / 2.8) +
    scale_color_manual(values = c("#CDCDCD", "#636ab1")) +
    scale_y_continuous(limits = c(0, 1), name = "Mean SMF Signal") +
    scale_x_continuous(name = "Relative position", limits = c(0, width),
      breaks = if (width == 512) c(0, 250, 500) else c(0, 512, 1024))
  path <- file.path(output_dir, filename)
  pdf(path, width = 20 / 2.54, height = 20 / 2.54)
  grid.arrange(set_panel_size(plot, height = unit(1.5, "cm"), width = unit(2, "cm")))
  dev.off()
  message("Wrote: ", path)
}

#Use the same more-than-five-TFBS criterion as the triangle plots.
for (region in unique(co_ac$region)) {
  if (nrow(motifs_for_region(region)) <= 5) next
  save_average(region, paste0(gsub("[:-]", "_", region), "_average.pdf"))
}

#Retain the original 512 bp Figure E3f example window without a GenomicRanges dependency.
example_region <- "chr1:132379834-132381882"
if (example_region %in% co_ac$region) {
  save_average(example_region, "Fig1_example_chr1_132380596_132381108_average.pdf",
               offset = 132380596 - (132379834 + 512), width = 512, motif_half_width = 2)
}
