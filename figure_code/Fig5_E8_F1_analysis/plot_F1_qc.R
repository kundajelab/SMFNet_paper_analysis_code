#Run from the repository root: Rscript figure_code/Fig5_E8_F1_analysis/plot_F1_qc.R
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
library(ggrastr)

PLOT_DIR <- "output/Fig5_E8_F1_analysis/plots/"

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6), 
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

read_lengths <- read.delim("input_files/Fig5_E8_F1_analysis/read_lengths_for_plotting.tsv")
n50_df <- read.delim("input_files/Fig5_E8_F1_analysis/read_length_N50.tsv") %>% 
  unite("sample", c(replicate, haplotype))

#Plot read lengths as violin
len_plot <- read_lengths %>% 
  unite("sample", c(replicate, haplotype)) %>% 
  ggplot() +
  geom_violin(aes(x = sample, y = read_length), fill = "#DEDEDE") +
  stat_summary(data = n50_df, aes(x = sample, y = N50), geom = "crossbar", lwd = 0.5, width = 0.25) +
  geom_text(data = n50_df, aes(x = sample, y = 40000, label = paste0("N50 = ", round(N50, 0))),
            inherit.aes = FALSE,size = 6/2.8) +
  scale_y_continuous(limits = c(0, 50000))

pdf(paste0(PLOT_DIR, "length_qc_f1.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(len_plot, height = unit(2.5, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
dev.off()

#Compute pair-wise correlation between haplotypes/reps
pos_fractions_reps <- fread(
  "input_files/Fig5_E8_F1_analysis/smf_per_position_per_sample.tsv",
  showProgress = TRUE)

#Calculate average coverage
cov_filter <- pos_fractions_reps %>% 
  select(sample, chr, start, end, coverage) %>% 
  filter(coverage >= 40) %>% 
  pivot_wider(names_from = sample, values_from = coverage) %>% 
  na.omit() %>% 
  select(chr, start, end)

fil_pos_frations <- pos_fractions_reps %>% 
  inner_join(cov_filter)

#The six unique comparisons in the orientation shown in Fig. E8c.
reps <- c("r1_BL6", "r1_SPRET", "r2_BL6", "r2_SPRET")
for (pair in combn(reps, 2, simplify = FALSE)) {
  rep_a <- pair[2]
  rep_b <- pair[1]
  loop_df <- fil_pos_frations %>% 
    filter(sample %in% c(rep_a, rep_b)) %>%
    select(-coverage, -smf_percent) %>% 
    pivot_wider(names_from = sample, values_from = smf_fraction)
  
  cor_loop <- cor(loop_df[[rep_a]], loop_df[[rep_b]])
  
  cor_plot <- loop_df %>% 
    ggplot(aes(x = .data[[rep_a]], y = .data[[rep_b]])) +
    rasterize(geom_hex(bins = 300), dpi = 1200) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_viridis(option = "A", trans = "log10", limits = c(1,100), oob = scales::squish) +
    annotate(geom = "text", label = paste0("r=", round(cor_loop, 2)), size = 5/2.8, x = 0.2, y = 0.9)
  
  
  pdf(paste0(PLOT_DIR, paste0("cor_reps_", rep_a, "_", rep_b, ".pdf")),
      width = 20/2.54, height = 20/2.54)
  fix <- set_panel_size(cor_plot, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
  grid.arrange(fix)
  dev.off()
}

#Compute correlation to ONT data
wt_comp <- fread("input_files/Fig5_E8_F1_analysis/smf_per_position_F1_WT.tsv")

fil_df <- wt_comp %>% 
  filter(F1_coverage >= 100 & WT_coverage >= 100)

cor_f1 <- cor(fil_df$F1_smf_fraction, fil_df$WT_smf_fraction)

cor_plot_wt <- fil_df %>% 
  ggplot(aes(x = F1_smf_fraction, y = WT_smf_fraction)) +
  rasterize(geom_hex(bins = 300), dpi = 1200) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_viridis(option = "A", trans = "log10") +
  annotate(geom = "text", label = paste0("r=", round(cor_f1, 2)), size = 5/2.8, x = 0.15, y = 0.95)

pdf(paste0(PLOT_DIR, paste0("cor_f1_wt_ont.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(cor_plot_wt, height = unit(3, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
dev.off()

