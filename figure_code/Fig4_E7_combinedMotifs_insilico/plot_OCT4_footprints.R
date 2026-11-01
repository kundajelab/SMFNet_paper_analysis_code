#Run from the repository root: Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_OCT4_footprints.R
library(data.table)
library(ggplot2)
library(egg)
library(gridExtra)
OUT <- 'output/Fig4_E7_combinedMotifs_insilico/oct4_footprints'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

#Read the matched BPNet and SMFNet summaries from compute_OCT4_footprints.py.
d <- fread(file.path(OUT, 'bpnet_delta_medians.tsv'))
r <- fread(file.path(OUT, 'bpnet_delta_correlations.tsv'))
s <- fread(file.path(OUT, 'threshold_summary.tsv'))
theme_set(theme_classic() + theme(
  panel.border = element_rect(color = 'black', fill = NA, linewidth = .5),
  axis.line = element_blank(),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 6),
  plot.title = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank(),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 6),
  plot.margin = margin(5, 5, 5, 5)))

#Add plotting margins around the observed changes.
pad <- function(v) {
  w <- diff(v)
  if (w == 0) w <- max(abs(v), 1) * .2
  v + c(-.18, .28) * w
}

#Compare model changes across spacings separately for KLF and ESRRB flanks.
plots <- list()
for (partner_name in c('SP-KLF', 'ESRRB')) {
  z <- d[partner == partner_name]
  rr <- r[partner == partner_name, pearson_r]
  p <- ggplot(z, aes(bpnet_median, smfnet_median)) +
    geom_point(aes(fill = spacing), shape = 21, color = 'black', size = 1.4, stroke = .3) +
    scale_fill_gradient(low = 'white', high = 'black', limits = c(20, 100), breaks = c(20, 60, 100),
      name = 'Spacing (bp)', guide = guide_colorbar(barheight = unit(1.5, 'cm'),
      barwidth = unit(.18, 'cm'))) +
    annotate('text', x = Inf, y = Inf, label = sprintf('Pearson r = %.3f', rr), hjust = 1.05,
      vjust = 1.3, size = 2) +
    coord_cartesian(xlim = pad(range(z$bpnet_median)), ylim = pad(range(z$smfnet_median)), expand = FALSE) +
    scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 3) +
    labs(x = 'BPNet count change', y = 'SMFNet bound change (pp)',
      title = paste(if (partner_name == 'SP-KLF') 'KLF' else 'ESRRB', '|', 'all molecules'))
  plots[[partner_name]] <- set_panel_size(p, width = unit(2.5, 'cm'), height = unit(2.5, 'cm'))
}
title <- grid::textGrob('Strict span <70 bp | change from center-only OCT4::SOX2',
  gp = grid::gpar(fontsize = 7))
pdf(file.path(OUT, 'BPNet_vs_SMFNet_delta.pdf'), width = 12 / 2.54, height = 7 / 2.54)
grid.arrange(grobs = plots, ncol = 2, top = title, newpage = FALSE)
dev.off()

#Plot the binding change with bootstrap intervals for all versus accessible molecules.
colors <- c(all_binding = '#646db3', open_binding = '#D55E00')
plots <- list()
for (partner_name in c('SP-KLF', 'ESRRB')) {
  z <- s[partner == partner_name]
  p <- ggplot(z, aes(spacing, value, color = metric, fill = metric)) +
    geom_hline(yintercept = 0, color = '#777777', linetype = 'dashed', linewidth = .3) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = .18,
      color = NA) + geom_line(linewidth = .35) + geom_point(size = .6) +
    scale_color_manual(values = colors, breaks = names(colors), labels = c('All molecules',
      'Accessible only')) +
    scale_fill_manual(values = colors, breaks = names(colors), labels = c('All molecules',
      'Accessible only')) +
    scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
    labs(x = 'Spacing (bp)', y = 'Bound fraction change (pp)',
      title = if (partner_name == 'SP-KLF') 'KLF flanks' else 'ESRRB flanks') +
    theme(legend.position = 'bottom', legend.title = element_blank()) +
    guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
  plots[[partner_name]] <- set_panel_size(p, width = unit(2.5, 'cm'), height = unit(1.25, 'cm'))
}
pdf(file.path(OUT, 'spacing_vs_delta_bound_fraction.pdf'), width = 12 / 2.54, height = 7 / 2.54)
grid.arrange(grobs = plots, ncol = 2, top = title, newpage = FALSE)
dev.off()
