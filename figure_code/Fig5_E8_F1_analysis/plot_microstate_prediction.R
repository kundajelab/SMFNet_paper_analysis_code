#Run from the repository root: Rscript figure_code/Fig5_E8_F1_analysis/plot_microstate_prediction.R
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
library(hdf5r)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6), 
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

PLOT_DIR <- "output/Fig5_E8_F1_analysis/plots/"
cluster_df <- read.delim("output/Fig5_E8_F1_analysis/F1_all4_balanced_cluster_fractions.tsv") %>% 
  select(-any_of(c("k", "silhouette", "count", "total_reads"))) %>% 
  pivot_wider(names_from = read_source, values_from = fraction, values_fill = 0) %>% 
  mutate(actual_delta = actual_b6 - actual_spret, predicted_delta = predicted_b6 - predicted_spret) 

out_table <- read.delim("output/Fig5_E8_F1_analysis/F1_all4_balanced_cluster_fractions.tsv") %>%
  select(-any_of(c("k", "silhouette", "count"))) %>%
  pivot_wider(
    names_from = read_source,
    values_from = c(fraction, total_reads),
    names_glue = "{read_source}_{.value}",
    values_fill = 0
  ) %>%
  mutate(actual_delta = actual_b6_fraction - actual_spret_fraction,
    predicted_delta = predicted_b6_fraction - predicted_spret_fraction)

write.table(out_table, "output/Fig5_E8_F1_analysis/f1_pred_obs_cluster_delta.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cluster_cor <- cor(cluster_df$actual_delta, cluster_df$predicted_delta)

cluster_delta <- cluster_df %>% 
  ggplot(aes(x = actual_delta, y = predicted_delta)) +
  geom_hex(bins = 100) +
  scale_fill_viridis(option = "A", trans = "log10") +
  coord_cartesian(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6)) +
  annotate(geom = "text", x = -0.45, y = 0.5, label = paste0("r=", round(cluster_cor, 2)), size = 5/2.8) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") 

pdf(paste0(PLOT_DIR,"delta_cluster_frac.pdf"),width=20/2.54,height=20/2.54)
grid.arrange(set_panel_size(cluster_delta,height=unit(3,"cm"),width=unit(3,"cm"))); dev.off()

