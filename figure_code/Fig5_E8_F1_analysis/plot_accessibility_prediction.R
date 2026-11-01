#Run from the repository root: Rscript figure_code/Fig5_E8_F1_analysis/plot_accessibility_prediction.R
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

ac_df <- read.delim("output/Fig5_E8_F1_analysis/F1_ac_pred_act.tsv")

fil_ac_df <- ac_df %>% 
  na.omit() %>% 
  filter(actual_b6_n >= 50) %>% 
  filter(actual_spret_n >= 50)

write.table(fil_ac_df, "output/Fig5_E8_F1_analysis/f1_pred_obs_ac_delta.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

delta_cor <- cor(fil_ac_df$actual_delta_accessible, fil_ac_df$predicted_delta_accessible)

delta_plot <- fil_ac_df %>% 
  ggplot() +
  geom_hex(aes(x = actual_delta_accessible, y = predicted_delta_accessible), bins = 60) +
  scale_fill_viridis(option = "A", transform = "log10") +
  annotate(geom = "text", x = -0.35, y = 0.45, label = paste0("r=", round(delta_cor, 2)), size = 5/2.8) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") 

pdf(paste0(PLOT_DIR,"delta_ac_perc_pred.pdf"),width=20/2.54,height=20/2.54)
grid.arrange(set_panel_size(delta_plot,height=unit(2,"cm"),width=unit(2,"cm"))); dev.off()

