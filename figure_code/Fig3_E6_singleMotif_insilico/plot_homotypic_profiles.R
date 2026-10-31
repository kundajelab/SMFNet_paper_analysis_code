#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig3_E6_singleMotif_insilico/plot_homotypic_profiles.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(egg)
library(gridExtra)
library(data.table)
library(viridis)


theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6),
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))
BIN_PLOT_DIR <- "output/Fig3_E6_singleMotif_insilico/homotypic_bins/"

dir.create(BIN_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

avg_sum_df <- read.delim("output/Fig3_E6_singleMotif_insilico/homotypic_spacing50_average_profiles.tsv") %>%
  separate(combo, c("spacing", "copies", "motif"), sep = "_") %>%
  filter(spacing == 50) %>%
  transmute(motif, copies = str_extract(copies, "[0-9]"), pos = position, signal = 1 - openness)

nomotif_sum <- avg_sum_df %>%
  filter(motif== "noMotif") %>%
  filter(between(pos, 128, 896)) %>%
  mutate(pos = pos - 512)

bin_df <- fread("output/Fig3_E6_singleMotif_insilico/homotypic_spacing50_binned_profiles.tsv")%>%
  separate(combo, c("spacing", "copies", "motif"), sep = "_") %>%
  filter(spacing == 50) %>%
  transmute(motif, copies = str_extract(copies, "[0-9]"), pos = position, signal = 1 - openness, bin = pseudo_bin + 1)

combos <- bin_df %>%
  distinct(motif, copies)

for (i in seq_len(nrow(combos))) {

  m <- combos$motif[i]
  cp <- combos$copies[i]

  print(paste0("Plotting ", m, " ", cp, " copies"))

  safe_m <- gsub("::", "_", m)
  safe_m <- gsub("-", "_", safe_m)

  sum_loop <- avg_sum_df %>%
    filter(motif == m, copies == cp) %>%
    filter(between(pos, 128, 896)) %>%
    mutate(pos = pos - 512)

  sum_plot <- sum_loop %>%
    ggplot() +
    geom_line(data = nomotif_sum, aes(x = pos, y = signal), color = "#BCBCBC") +
    geom_line(aes(x = pos, y = signal), color = "#646db3") +
    scale_x_continuous(expand = c(0, 0), breaks = c(-256, 0, 256)) +
    scale_y_continuous(limits = c(0, 1))

  pdf(paste0(BIN_PLOT_DIR, "", safe_m, "_", cp, "homotypic_motif_avg.pdf"),
      width = 20/2.54, height = 20/2.54)
  fix <- set_panel_size(sum_plot, height = grid::unit(1.5, "cm"), width = grid::unit(1.5, "cm"))
  grid.arrange(fix, newpage = FALSE)
  dev.off()


  bin_loop <- bin_df %>%
    filter(motif == m, copies == cp) %>%
    filter(between(pos, 128, 896)) %>%
    mutate(pos = pos - 512)

  bin_plot <- bin_loop %>%
    ggplot() +
    geom_raster(aes(x = pos, y = factor(bin), fill = signal)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(-256, 0, 256)) +
    scale_y_discrete(expand = c(0, 0), name = "Accessibility bin", breaks = c(1, 5, 10, 15)) +
    scale_fill_viridis(option = "A", limits = c(0, 1), direction = -1) +
    theme(axis.ticks.y = element_blank())

  pdf(paste0(BIN_PLOT_DIR, "", safe_m, "_", cp, "homotypic_motif_bin.pdf"),
      width = 20/2.54, height = 20/2.54)
  fix <- set_panel_size(bin_plot, height = grid::unit(2, "cm"), width = grid::unit(1.5, "cm"))
  grid.arrange(fix, newpage = FALSE)
  dev.off()
}

