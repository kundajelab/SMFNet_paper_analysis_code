#Run from the repository root: Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_state_profiles.R
library(data.table)
library(ggplot2)
library(egg)
library(gridExtra)
library(viridis)
OUT <- "output/Fig4_E7_combinedMotifs_insilico/state_profiles"

#Read pooled average and binned profiles from compute_state_profiles.py.
avg <- fread(file.path(OUT, "average_profiles.tsv"))
bins <- fread(file.path(OUT, "binned_profiles.tsv"))

#Crop to the central window and express positions relative to the motif center.
avg <- avg[position >= 128 & position <= 896][, `:=`(pos = position-512, signal = 1-openness)]
bins <- bins[position >= 128 & position <= 896][, `:=`(pos = position-512, signal = 1-openness,
  bin = pseudo_bin + 1)]
theme_set(theme_classic() + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
  axis.line = element_blank(),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 6)))
save_panel <- function(p, name, height) {
  pdf(file.path(OUT, paste0(name, ".pdf")), width = 20 / 2.54, height = 20 / 2.54)
  grid.arrange(set_panel_size(p, width = unit(1.5, "cm"), height = unit(height, "cm")), newpage = FALSE)
  dev.off()
}

#Write an average profile and an accessibility-bin heatmap for each selected construct.
combos <- unique(bins[, .(analysis, spacing, combo)])
combos <- combos[analysis == "enhancer" | !grepl("noMotif_OCT4::SOX2_noMotif$", combo)]
for (i in seq_len(nrow(combos))) {
  analysis_name <- combos$analysis[i]
  s <- combos$spacing[i]
  combination <- combos$combo[i]
  z <- avg[analysis == analysis_name & spacing == s & combo == combination]
  #Use the original center-only profile used as the grey overlay in the paper.
  control <- avg[analysis == analysis_name & spacing == s & combo == paste0(s, "_noMotif_OCT4::SOX2_noMotif")]
  p <- ggplot() + geom_line(data = control, aes(pos, signal), color = "#BCBCBC") +
    geom_line(data = z, aes(pos, signal), color = "#646db3") +
    scale_x_continuous(expand = c(0, 0), breaks = c(-256, 0, 256)) +
    scale_y_continuous(limits = c(0, 1))
  name <- paste0(analysis_name, "_", gsub("::|-", "_", combination))
  save_panel(p, paste0(name, "_average"), 1.5)
  #Display the same state-3 signal in each read bin as in the average profiles.
  z <- bins[analysis == analysis_name & spacing == s & combo == combination]
  p <- ggplot(z) + geom_raster(aes(pos, factor(bin), fill = signal)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(-256, 0, 256)) +
    scale_y_discrete(expand = c(0, 0), name = "Accessibility bin", breaks = c(1, 5, 10, 15)) +
    scale_fill_viridis(option = "A", limits = c(0, 1), direction = -1) + theme(axis.ticks.y = element_blank())
  save_panel(p, paste0(name, "_bins"), 2)
}
