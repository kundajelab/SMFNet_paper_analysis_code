#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig3_E6_singleMotif_insilico/plot_single_motifs.R

INPUT_DIR <- "input_files/Fig3_E6_singleMotif_insilico"
INTERMEDIATE_DIR <- "output"
RESULTS_DIR <- "output"
library(dplyr)
library(tidyr)
library(ggplot2)
library(egg)
library(gridExtra)
library(viridis)
library(ggh4x)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6),
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

PLOT_DIR <- paste0(file.path(RESULTS_DIR, "Fig3_E6_singleMotif_insilico"), "/")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
df <- read.delim(file.path(INTERMEDIATE_DIR, "Fig3_E6_singleMotif_insilico/single_motif_spacing_50bp_ac_summary_median_dispersion_noMin_width.txt"))

fil_df <- df %>%
  separate(region, c("bg_region", "spacing", "motif1", "motif2", "motif3"), "_") %>%
  filter(motif1 == "noMotif" & motif3 == "noMotif") %>%
  transmute(motif = motif2, bg_region, width = med_width, freq = perc_ac, dispersion = boundary_dispersion)

cor(fil_df$freq, fil_df$dispersion)
cor(fil_df$width, fil_df$dispersion)
cor(fil_df$width, fil_df$freq)

plot_df <- fil_df %>%
  pivot_longer(-c(motif, bg_region), names_to = "measure", values_to = "value")

line_df <- plot_df %>%
  filter(motif == "noMotif") %>%
  group_by(measure) %>%
  summarize(median_val = median(value))

motif_order <- plot_df %>%
  filter(measure == "freq") %>%
  group_by(motif) %>%
  summarize(frac = median(value)) %>%
  arrange(if_else(motif == "noMotif", 0, 1), frac) %>%
  pull(motif)

gradient_col <- colorRampPalette(c("#ece2f0", "#68a9cf", "#056551"))(20)

sum_df <- plot_df %>%
  group_by(measure, motif) %>%
  summarise(
    med = median(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    n_regions = n(),
    .groups = "drop"
  )

write.table(sum_df, file.path(INTERMEDIATE_DIR, "Fig3_E6_singleMotif_insilico/single_motif_summary.tsv"),
            sep = "\t", quote = FALSE)

simple_bars <- sum_df %>%
  ggplot(aes(y = factor(motif, levels = rev(motif_order)), fill = factor(motif, levels = motif_order))) +
  facet_wrap(~measure, scales = "free", nrow = 1) +
  geom_errorbarh(aes(xmin = q1, xmax = q3), height = 0, linewidth = 0.7) +
  geom_point(aes(x = med), shape = 21, size = 2, color = "black") +
  facetted_pos_scales(x = list(
    measure == "width" ~ scale_x_continuous(limits = c(110, 260), breaks = c(150, 200, 250)),
    measure == "dispersion" ~ scale_x_continuous(limits = c(70, 200), breaks = c(100, 150, 200)),
    measure == "freq" ~ scale_x_continuous(limits = c(0.2, 0.9), breaks = c(0.25, 0.5, 0.75))
  )) +
  scale_fill_manual(values = c("#BCBCBC", gradient_col)) +
  geom_vline(data = line_df, aes(xintercept = median_val), linetype = "dashed") +
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.title = element_blank())

pdf(paste0(PLOT_DIR, paste0("singleMotif_measures.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(simple_bars, height = grid::unit(4.5, "cm"), width = grid::unit(2, "cm"))
grid.arrange(fix)
dev.off()

dis_freq_cor <- cor(sum_df[sum_df$measure == 'dispersion',]$med, sum_df[sum_df$measure == 'freq',]$med)
dis_wid_cor <- cor(sum_df[sum_df$measure == 'dispersion',]$med, sum_df[sum_df$measure == 'width',]$med)
wid_freq_cor <- cor(sum_df[sum_df$measure == 'width',]$med, sum_df[sum_df$measure == 'freq',]$med)

dis_freq_scatter <- sum_df %>%
  select(measure, motif, med) %>%
  pivot_wider(names_from = measure, values_from = med) %>%
  ggplot() +
  geom_point(aes(x = dispersion, y = freq, fill = factor(motif, levels = motif_order)), shape = 21) +
  geom_text(aes(x = dispersion, y = freq, label = motif), size = 5/2.8) +
  scale_fill_manual(values = c("#BCBCBC", gradient_col)) +
  annotate(geom = "text", x = 130, y = 0.8, label = paste0("r=", round(dis_freq_cor, 2)), size = 6/2.8) +
  theme(legend.position = "none") +
  geom_smooth(aes(x = dispersion, y = freq), method = "lm", color = "black", se = FALSE, linetype = "dashed") +
  scale_x_continuous(limits = c(80, 160)) +
  scale_y_continuous(limits = c(0.25, 0.85))

pdf(paste0(PLOT_DIR, paste0("singleMotif_dis_freq_scatter.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(dis_freq_scatter, height = grid::unit(3, "cm"), width = grid::unit(3, "cm"))
grid.arrange(fix)
dev.off()

BIN_PLOT_DIR <- paste0(file.path(RESULTS_DIR, "Fig3_E6_singleMotif_insilico"), "/")

sum_df <- read.delim(file.path(INTERMEDIATE_DIR, "Fig3_E6_singleMotif_insilico/example_three_motif_spacing_50bp_average_profiles.tsv")) %>%
  separate(combo, c("spacing", "motif1", "motif2", "motif3"), sep = "_") %>%
  filter(motif1 == "noMotif" & motif3 == "noMotif") %>%
  transmute(motif = motif2, pos = position + 1, signal = 1 - openness)

nomotif_sum <- sum_df %>%
  filter(motif== "noMotif") %>%
  filter(between(pos, 128, 896)) %>%
  mutate(pos = pos - 512)

bin_df <- read.delim(file.path(INTERMEDIATE_DIR, "Fig3_E6_singleMotif_insilico/example_three_motif_spacing_50bp_binned_profiles.tsv"))%>%
  separate(combo, c("spacing", "motif1", "motif2", "motif3"), sep = "_") %>%
  filter(motif1 == "noMotif" & motif3 == "noMotif") %>%
  transmute(motif = motif2, pos = position + 1, signal = 1 - openness, bin = pseudo_bin + 1)

motifs <- unique(bin_df$motif)

for (m in motifs) {
  print(paste0("Plotting ", m))

  safe_m <- gsub("::", "_", m)
  safe_m <- gsub("-", "_", safe_m)

  sum_loop <- sum_df %>%
    filter(motif == m) %>%
    filter(between(pos, 128, 896)) %>%
    mutate(pos = pos - 512)

  sum_plot <- sum_loop %>%
    ggplot() +
    geom_line(data = nomotif_sum, aes(x = pos, y = signal), color = "#BCBCBC") +
    geom_line(aes(x = pos, y = signal), color = "#646db3") +
    scale_x_continuous(expand = c(0, 0), breaks = c(-256, 0, 256)) +
    scale_y_continuous(limits = c(0,1))

  pdf(paste0(BIN_PLOT_DIR, "", safe_m, "_solo_motif_avg.pdf"),
      width = 20/2.54, height = 20/2.54)
  fix <- set_panel_size(sum_plot, height = grid::unit(1.5, "cm"), width = grid::unit(1.5, "cm"))
  grid.arrange(fix)
  dev.off()

  bin_loop <-  bin_df %>%
    filter(motif == m) %>%
    filter(between(pos, 128, 896)) %>%
    mutate(pos = pos - 512)

  bin_plot <- bin_loop %>%
    ggplot() +
    geom_raster(aes(x = pos, y = factor(bin), fill = signal)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(-256, 0, 256)) +
    scale_y_discrete(expand = c(0, 0), name = "Accessibility bin", breaks = c(1,5,10,15)) +
    scale_fill_viridis(option = "A", limits = c(0, 1), direction = -1) +
    theme(axis.ticks.y = element_blank())

  pdf(paste0(BIN_PLOT_DIR, "", safe_m, "_solo_motif_bin.pdf"),
      width = 20/2.54, height = 20/2.54)
  fix <- set_panel_size(bin_plot, height = grid::unit(2, "cm"), width = grid::unit(1.5, "cm"))
  grid.arrange(fix)
  dev.off()

}
