#Run from the repository root: Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_phasing.R
library(data.table)
library(ggplot2)
library(egg)
library(gridExtra)
OUT <- 'output/Fig4_E7_combinedMotifs_insilico/median_global_residuals'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

#Read matched-background changes, global fit coefficients and bootstrap residual intervals.
d <- fread(file.path(OUT, 'all_spacings_scatter_data.txt'))
r <- fread(file.path(OUT, 'OCT4_SOX2_global_residual_spacing_CIs.txt'))
cf <- fread(file.path(OUT, 'global_all_spacings_linear_model.txt'))
a <- cf[term == 'intercept', estimate]
b <- cf[term == 'slope', estimate]
stopifnot(max(abs(d$residual-(d$dispersion_delta-a-b * d$freq_delta)))<1e-10)
COLOR <- '#D55E00'
GREY <- '#666666'
focus <- c('SP-KLF', 'ESRRB', 'OCT4::SOX2')
r[, flank_motif := factor(flank_motif, levels = focus)]
theme_set(theme_classic() + theme(
  panel.border = element_rect(color = 'black', fill = NA, linewidth = .5),
  axis.line = element_blank(),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 6),
  title = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank(),
  plot.caption = element_text(size = 6)))
save_plot <- function(p, name, w, h) {
  fix <- set_panel_size(p, width = unit(w, 'cm'), height = unit(h, 'cm'))
  pdf(file.path(OUT, paste0(name, '.pdf')), width = 20 / 2.54, height = 20 / 2.54)
  grid.arrange(fix, newpage = FALSE)
  dev.off()
}

#Compare median frequency and dispersion changes across all retained motif combinations.
pearson_all <- cor(d$freq_delta, d$dispersion_delta, method = 'pearson')
fwrite(data.table(scope = 'all', n = nrow(d), pearson_r = pearson_all), file.path(OUT,
  'scatter_pearson.txt'), sep = '\t')
p <- ggplot(d, aes(freq_delta, dispersion_delta)) +
  annotate('text', x = Inf, y = Inf, label = sprintf('Pearson r = %.3f', pearson_all), hjust = 1.08,
    vjust = 1.4, size = 2) +
  geom_hline(yintercept = 0, color = 'grey85', linewidth = .2) + geom_vline(xintercept = 0,
    color = 'grey85', linewidth = .2) +
  geom_point(color = GREY, size = .45, alpha = .35, stroke = 0) +
  geom_abline(intercept = a, slope = b, color = COLOR, linetype = 'dashed', linewidth = .5) +
  labs(x = 'Median frequency delta', y = 'Median dispersion delta (bp)',
    title = 'All retained combinations | 20-100 bp',
    caption = 'Matched-background deltas, then medians; dashed orange: global linear fit')
save_plot(p, 'dispersion_vs_frequency_all_spacings_grey', 3, 3)

#Plot residual dispersion after subtracting the global frequency relationship.
p <- ggplot(r, aes(spacing, residual)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = '#646db3', alpha = .2) +
  geom_line(color = '#646db3', linewidth = .35) + geom_point(color = '#646db3', size = .6) +
  geom_hline(yintercept = 0, color = COLOR, linetype = 'dashed', linewidth = .5) +
  facet_wrap(~flank_motif, ncol = 1, scales = 'fixed') +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
  labs(x = 'Spacing (bp)', y = 'Dispersion residual (bp)', title = 'OCT4::SOX2 centre',
    caption = '95% background-bootstrap CI; medians and global fit recomputed each draw\nDashed orange zero corresponds to the global regression line in the scatter')
save_plot(p, 'OCT4_SOX2_global_residual_vs_spacing', 2.5, 1.25)

#Boundary-dispersion changes by spacing (E7k).
p <- ggplot(d, aes(spacing, dispersion_delta, group = spacing, fill = spacing)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3, color = "#777777") +
  geom_boxplot(outlier.shape = NA, width = 3.5, linewidth = .3, color = "#245780") +
  scale_fill_gradient(low = "white", high = "black", limits = c(20, 100)) +
  scale_x_continuous(breaks = seq(20, 100, 20), limits = c(16, 104), expand = c(0, 0)) +
  labs(x = "Motif spacing (bp)", y = "Boundary dispersion change (bp)") + theme(legend.position = "none")
save_plot(p, "phasing_by_spacing", 3.5, 2)
