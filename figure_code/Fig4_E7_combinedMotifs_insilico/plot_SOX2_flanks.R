#Run from the repository root: Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_SOX2_flanks.R
library(ggplot2)
library(dplyr)
library(tibble)
library(egg)
library(gridExtra)
library(data.table)
library(viridis)
filter <- dplyr::filter
summarize <- dplyr::summarize

OUT <- 'output/Fig4_E7_combinedMotifs_insilico/SOX2_flanks'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
INPUT <- 'output/Fig4_E7_combinedMotifs_insilico/selected_pairs_baseline'
setDTthreads(4)
theme_set(theme_classic() +
  theme(
  legend.text = element_text(size = 6),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(),
  axis.text.x = element_text(size = 6),
  title = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.title = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank()))

#Compare SOX2-flanked constructs with their matched center-only backgrounds.
quant <- fread(file.path(INPUT, "combined_quantification.txt"),
  select = c("spacing", "bg_region", "left", "center", "right", "perc_ac"))
sox2 <- quant[left == "SOX2" & right == "SOX2" & center %in% c("OCT4::SOX2", "SOX2", "ESRRB", "SP-KLF")]
center_only <- quant[left == "noMotif" & right == "noMotif",
.(spacing, bg_region, center, center_frequency = perc_ac)]
paired <- merge(sox2, center_only, by = c("spacing", "bg_region", "center"), all.x = TRUE)

#Average the paired frequency change across backgrounds for each center motif and spacing.
delta_freq_df <- as_tibble(paired) %>%
  mutate(delta_freq = perc_ac-center_frequency) %>%
  group_by(center_motif = center, spacing) %>%
  summarize(delta_freq = mean(delta_freq), .groups = "drop")
fwrite(delta_freq_df, file.path(OUT, "SOX2_flanks_frequency.tsv"), sep = "\t")
PLOT_DIR <- paste0(OUT, "/")

#Show the spacing dependence with the original smoothed trend.
delta_freq_plot <- delta_freq_df %>%
  filter(center_motif %in% c("OCT4::SOX2", "SOX2", "ESRRB", "SP-KLF")) %>%
  ggplot(aes(x = as.numeric(spacing), y = delta_freq)) +
  facet_wrap(~center_motif) +
  geom_point(color = "#646db3") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_smooth(se = FALSE, color = "black") +
  scale_y_continuous(limits = c(-0.05, 0.1))

pdf(paste0(PLOT_DIR, paste0("SOX2_flanks_frequency_LOESS.pdf")),
  width = 20 / 2.54, height = 20 / 2.54)
fix <- set_panel_size(delta_freq_plot, height = unit(1.5, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
dev.off()
