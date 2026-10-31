#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig3_E6_singleMotif_insilico/plot_gc_motifs.R

INPUT_DIR <- "input_files/Fig3_E6_singleMotif_insilico"
INTERMEDIATE_DIR <- "output"
RESULTS_DIR <- "output"
PLOT_DIR <- paste0(file.path(RESULTS_DIR, "Fig3_E6_singleMotif_insilico"), "/")
dir.create(PLOT_DIR, recursive=TRUE, showWarnings=FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(egg)
library(gridExtra)
library(viridis)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6),
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

df <- read.delim(file.path(INTERMEDIATE_DIR, "Fig3_E6_singleMotif_insilico/gc_single_motif_ac_summary.txt")) %>%
  separate(region, c("bg_region", "GC", "trash", "motif"), sep = "_") %>%
  transmute(motif, bg_region, GC = str_extract(GC, "[0-9][0-9]"), width = med_width, freq = perc_ac,
            dispersion = boundary_dispersion)

med_df <- df %>%
  group_by(motif, GC) %>%
  summarize(width = median(width), freq = median(freq), dispersion = mean(dispersion), .groups = "drop")

mat <- med_df %>%
  filter(motif != "noMotif") %>%
  select(motif, GC, freq) %>%
  tidyr::pivot_wider(names_from = GC, values_from = freq) %>%
  tibble::column_to_rownames("motif") %>%
  as.matrix()

freq_order <- rownames(mat)[order(rowMeans(mat))]

freq_heat <- med_df %>%
  ggplot() +
  geom_tile(aes(x = GC, y = factor(motif, levels = rev(c("noMotif", freq_order))), fill = freq)) +
  scale_fill_viridis(option = "A") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(breaks = c(30, 35, 40, 45, 50, 55, 60, 65, 70), expand = c(0,0), name = "Background %GC") +
  theme(axis.ticks = element_blank(), axis.title.y = element_blank())

pdf(paste0(PLOT_DIR, paste0("gc_freq_heat.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(freq_heat, height = grid::unit(4, "cm"), width = grid::unit(2.5, "cm"))
grid.arrange(fix)
dev.off()

med_df_z <- med_df %>%
  group_by(motif) %>%
  mutate(
    freq_z = (freq - mean(freq, na.rm = TRUE)) / sd(freq, na.rm = TRUE)
  ) %>%
  ungroup()

freq_z_heat <- med_df_z %>%
  ggplot() +
  geom_tile(aes(x = GC, y = factor(motif, levels = rev(c("noMotif", freq_order))), fill = freq_z)) +
  scale_fill_gradient2(low = "skyblue", high = "darkred", mid = "white", midpoint = 0, limits = c(-3, 3)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(breaks = c(30, 35, 40, 45, 50, 55, 60, 65, 70), expand = c(0,0), name = "Background %GC") +
  theme(axis.ticks = element_blank(), axis.title.y = element_blank())

pdf(paste0(PLOT_DIR, paste0("gc_freq_z_heat.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(freq_z_heat, height = grid::unit(4, "cm"), width = grid::unit(2.5, "cm"))
grid.arrange(fix)
dev.off()
