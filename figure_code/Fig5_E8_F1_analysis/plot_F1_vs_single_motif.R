#Run from the repository root: Rscript figure_code/Fig5_E8_F1_analysis/plot_F1_vs_single_motif.R
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
library(viridis)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6), 
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

PLOT_DIR <- "output/Fig5_E8_F1_analysis/plots/"

f1_metrics <- read.delim("output/Fig5_E8_F1_analysis/f1_metrics_sum.txt") %>% 
  transmute(motif, f1_width = med_width, f1_perc = med_perc, f1_disp = med_disp)

insilico_metrics <- read.delim("input_files/Fig5_E8_F1_analysis/insilico_metrics_sum.txt") %>%
  select(-q1, -q3, -n_regions) %>% 
  pivot_wider(names_from = measure, values_from = med) %>% 
  transmute(motif, silico_width = width, silico_perc = freq, silico_disp = dispersion)

combi_metrics <- inner_join(f1_metrics, insilico_metrics)

write.table(combi_metrics, 
            "output/Fig5_E8_F1_analysis/F1_vs_insilico_metrics.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

perc_cor <- cor(combi_metrics$f1_perc, combi_metrics$silico_perc)
width_cor <- cor(combi_metrics$f1_width, combi_metrics$silico_width)
disp_cor <- cor(combi_metrics$f1_disp, combi_metrics$silico_disp)

perc_cor_plot <- combi_metrics %>% 
  ggplot(aes(x = f1_perc, y = silico_perc)) +
  geom_point(alpha = 0.8, size = 0.8, color = "#676767") +
  scale_y_continuous(limits = c(0.2, 0.85)) +
  scale_x_continuous(limits = c(-0.06, 0.02)) +
  annotate(geom = "text", x = 0.005, y = 0.8, label = paste0("r=", round(perc_cor, 2)), size = 5/2.8) +
  geom_text(aes(label = motif), size = 5/2.8)
  

pdf(paste0(PLOT_DIR,"f1_metrics_perc_cor.pdf"), width=20/2.54, height=20/2.54)
grid.arrange(set_panel_size(perc_cor_plot, height=unit(3,"cm"), width=unit(3,"cm")))
dev.off()

width_cor_plot <- combi_metrics %>% 
  ggplot(aes(x = f1_width, y = silico_width)) +
  geom_point(alpha = 0.8, size = 0.8, color = "#676767") +
  scale_y_continuous(limits = c(140, 230)) +
  scale_x_continuous(limits = c(-40, 5)) +
  annotate(geom = "text", x = 0, y = 220, label = paste0("r=", round(width_cor, 2)), size = 5/2.8)+
  geom_text(aes(label = motif), size = 5/2.8)

pdf(paste0(PLOT_DIR,"f1_metrics_width_cor.pdf"), width=20/2.54, height=20/2.54)
grid.arrange(set_panel_size(width_cor_plot, height=unit(3,"cm"), width=unit(3,"cm")))
dev.off()

disp_cor_plot <- combi_metrics %>% 
  ggplot(aes(x = f1_disp, y = silico_disp)) +
  geom_point(alpha = 0.8, size = 0.8, color = "#676767") +
  scale_y_continuous(limits = c(80, 160)) +
  scale_x_continuous(limits = c(-30, 15)) +
  annotate(geom = "text", x = 8, y = 158, label = paste0("r=", round(disp_cor, 2)), size = 5/2.8) +
  geom_text(aes(label = motif), size = 5/2.8)

pdf(paste0(PLOT_DIR,"f1_metrics_disp_cor.pdf"), width=20/2.54, height=20/2.54)
grid.arrange(set_panel_size(disp_cor_plot, height=unit(3,"cm"), width=unit(3,"cm")))
dev.off()

