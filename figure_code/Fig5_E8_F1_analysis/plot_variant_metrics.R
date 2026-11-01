#Run from the repository root: Rscript figure_code/Fig5_E8_F1_analysis/plot_variant_metrics.R
dir.create("output/Fig5_E8_F1_analysis/plots", recursive=TRUE, showWarnings=FALSE)
PLOT_DIR <- "output/Fig5_E8_F1_analysis/plots/"

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(forcats)
filter <- dplyr::filter
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

PLOT_DIR <- "output/Fig5_E8_F1_analysis/plots/"

df <- read.delim("output/Fig5_E8_F1_analysis/F1_atacPeaks_variant_hits_ac_metrics_median_disp.tsv")

fil_df <- df %>% 
  mutate(cwm_delta = source_cwm_score - opposite_cwm_score) %>% 
  mutate(delta_perc_ac = ifelse(source_allele=="SPRET", -delta_perc_ac, delta_perc_ac),
         delta_median_width = ifelse(source_allele=="SPRET", -delta_median_width, delta_median_width),
         delta_boundary_dispersion = ifelse(source_allele=="SPRET", -delta_boundary_dispersion, delta_boundary_dispersion)) %>% 
  group_by(motif) %>% 
  slice_max(order_by = cwm_delta, prop = 0.1, with_ties = FALSE) %>%
  ungroup() %>% 
  group_by(region_id) %>% 
  slice_max(cwm_delta) %>%
  ungroup()

write.table(fil_df, 
            "output/Fig5_E8_F1_analysis/F1_atacPeaks_variant_hits_ac_metrics_filter.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#Filter for motifs that come somewhat frequent
n <- fil_df %>% count(motif)
fil_motifs <- n %>% filter(n>=50) %>%
  filter(!str_detect(motif, "alt")) %>%  filter(!str_detect(motif, "double")) %>% pull(motif)

sum_plot_df <- fil_df %>%
  filter(motif %in% fil_motifs) %>%
  group_by(motif) %>%
  summarise(
    med_width = median(delta_median_width, na.rm = TRUE),
    q1_width = quantile(delta_median_width, 0.25, na.rm = TRUE),
    q3_width = quantile(delta_median_width, 0.75, na.rm = TRUE),
    med_perc = median(delta_perc_ac, na.rm = TRUE),
    q1_perc = quantile(delta_perc_ac, 0.25, na.rm = TRUE),
    q3_perc = quantile(delta_perc_ac, 0.75, na.rm = TRUE),
    med_disp = median(delta_boundary_dispersion, na.rm = TRUE),
    q1_disp = quantile(delta_boundary_dispersion, 0.25, na.rm = TRUE),
    q3_disp = quantile(delta_boundary_dispersion, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(motif = fct_reorder(motif, -med_perc))

write.table(sum_plot_df, "output/Fig5_E8_F1_analysis/f1_metrics_sum.txt",
            sep = "\t", quote = FALSE)

motif_cols <- rev(grey.colors(length(fil_motifs), start = 0.05, end = 0.95))

plot_width <- sum_plot_df %>%
  ggplot(aes(y = motif, fill = motif)) +
  geom_errorbarh(aes(xmin = q1_width, xmax = q3_width), height = 0, linewidth = 0.7) +
  geom_point(aes(x = med_width), shape = 21, size = 2.5, color = "black") +
  scale_fill_manual(values = motif_cols) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-100, 40)) +
  theme(legend.position = "none", axis.ticks.y = element_blank())

plot_perc <- sum_plot_df %>%
  ggplot(aes(y = motif, fill = motif)) +
  geom_errorbarh(aes(xmin = q1_perc, xmax = q3_perc), height = 0, linewidth = 0.7) +
  geom_point(aes(x = med_perc), shape = 21, size = 2.5, color = "black") +
  scale_fill_manual(values = motif_cols) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-0.15, 0.1)) +
  theme(legend.position = "none", axis.ticks.y = element_blank())

plot_disp <- sum_plot_df %>%
  ggplot(aes(y = motif, fill = motif)) +
  geom_errorbarh(aes(xmin = q1_disp, xmax = q3_disp), height = 0, linewidth = 0.7) +
  geom_point(aes(x = med_disp), shape = 21, size = 2.5, color = "black") +
  scale_fill_manual(values = motif_cols) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", axis.ticks.y = element_blank())

pdf(paste0(PLOT_DIR,"delta_width_box.pdf"),width=20/2.54,height=20/2.54)
grid.arrange(set_panel_size(plot_width,height=unit(4,"cm"),width=unit(2.5,"cm"))); dev.off()

pdf(paste0(PLOT_DIR,"delta_perc_box.pdf"),width=20/2.54,height=20/2.54)
grid.arrange(set_panel_size(plot_perc,height=unit(4,"cm"),width=unit(2.5,"cm"))); dev.off()

pdf(paste0(PLOT_DIR,"delta_dispersion_box_median.pdf"),width=20/2.54,height=20/2.54)
grid.arrange(set_panel_size(plot_disp,height=unit(4,"cm"),width=unit(2.5,"cm"))); dev.off()

