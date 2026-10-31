#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig3_E6_singleMotif_insilico/plot_homotypic_motifs.R

INPUT_DIR <- "input_files/Fig3_E6_singleMotif_insilico"
INTERMEDIATE_DIR <- "output"
RESULTS_DIR <- "output"
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
PLOT_DIR <- paste0(file.path(RESULTS_DIR, "Fig3_E6_singleMotif_insilico"), "/")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
df <- read.delim(file.path(INTERMEDIATE_DIR, "Fig3_E6_singleMotif_insilico/homotypic_all_spacings_ac_summary_median_dispersion.txt"))

fil_df <- df %>%
  separate(region, c("bg_region", "spacing", "copies", "motif"), "_") %>%
  transmute(motif, bg_region, copies = as.numeric(str_extract(copies, "[0-9]")), spacing,
            width = med_width, freq = perc_ac, dispersion = boundary_dispersion)

sum_df <- fil_df %>%
  filter(spacing == 50) %>%
  group_by(motif, copies) %>%
  summarize(
    width_med = median(width, na.rm=TRUE),
    width_q1 = quantile(width, 0.25, na.rm=TRUE),
    width_q3 = quantile(width, 0.75, na.rm=TRUE),
    freq_med = median(freq, na.rm=TRUE),
    freq_q1 = quantile(freq, 0.25, na.rm=TRUE),
    freq_q3 = quantile(freq, 0.75, na.rm=TRUE),
    .groups="drop"
  )

motifs <- unique(sum_df$motif)

noMotif_df <- sum_df %>%
  filter(motif == "noMotif") %>%
  select(-motif) %>%
  tidyr::crossing(motif = motifs) %>%
  filter(motif != "noMotif")

wid_freq_plot <- rbind(noMotif_df, sum_df) %>%
  filter(motif %in% c("SP-KLF", "OCT4::SOX2", "CTCF", "ESRRB")) %>%
  ggplot(aes(x = width_med, y = freq_med, fill = as.factor(copies))) +
  facet_wrap(~motif) +
  geom_errorbar(aes(ymin=freq_q1, ymax=freq_q3), width=0) +
  geom_errorbarh(aes(xmin=width_q1, xmax=width_q3), height=0) +
  geom_point(size=2, shape = 21) +
  labs(x="Accessible width (bp)", y="Accessible fraction (%)") +
  scale_fill_grey(start = 1, end = 0, name = "Copies") +
  scale_y_continuous(name = "Accessible fraction (%)", limits = c(0, 1)) +
  scale_x_continuous(name = "Width (bp)", limits = c(100, 450))

pdf(paste0(PLOT_DIR, paste0("homotypic_width_freq.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(wid_freq_plot, height = grid::unit(1.5, "cm"), width = grid::unit(1.5, "cm"))
grid.arrange(fix)
dev.off()

med_df <- fil_df %>%
  filter(spacing == 50 & copies > 0) %>%
  pivot_longer(c(width, freq, dispersion), names_to = "measure", values_to = "val") %>%
  group_by(motif, copies, measure) %>%
  summarize(med = median(val, na.rm =TRUE),
            .groups="drop")

freq_heat <- med_df %>%
  filter(measure == "freq") %>%
  group_by(motif) %>%
  mutate(avg = mean(med, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(motif = reorder(motif, -avg)) %>%
  ggplot(aes(x=as.factor(copies), y=motif, fill=med)) +
  geom_tile(color = "#454545") +
  theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(expand = c(0,0), name = "Copy number") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_viridis(option = "A", limits = c(0.2,0.95))

pdf(paste0(PLOT_DIR, paste0("homotypic_freq_heat_50bp_spacing.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(freq_heat, height = grid::unit(4, "cm"), width = grid::unit(1.5, "cm"))
grid.arrange(fix)
dev.off()

width_heat <- med_df %>%
  filter(measure == "width") %>%
  group_by(motif) %>%
  mutate(avg = mean(med, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(motif = reorder(motif, -avg)) %>%
  ggplot(aes(x=as.factor(copies), y=motif, fill=med)) +
  geom_tile(color = "#454545") +
  theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(expand = c(0,0), name = "Copy number") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_viridis(option = "A", limits = c(140, 460))

pdf(paste0(PLOT_DIR, paste0("homotypic_width_heat_50bp_spacing.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(width_heat, height = grid::unit(4, "cm"), width = grid::unit(1.5, "cm"))
grid.arrange(fix)
dev.off()

dispersion_heat <- med_df %>%
  filter(measure == "dispersion") %>%
  group_by(motif) %>%
  mutate(avg = mean(med, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(motif = reorder(motif, avg)) %>%
  ggplot(aes(x=as.factor(copies), y=motif, fill=med)) +
  geom_tile(color = "#454545") +
  theme(axis.ticks = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(expand = c(0,0), name = "Copy number") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_viridis(option = "A", limits = c(70, 200), direction = -1, breaks = c(100, 150, 200))

pdf(paste0(PLOT_DIR, paste0("homotypic_dispersion_heat_50bp_spacing.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(dispersion_heat, height = grid::unit(4, "cm"), width = grid::unit(1.5, "cm"))
grid.arrange(fix)
dev.off()

