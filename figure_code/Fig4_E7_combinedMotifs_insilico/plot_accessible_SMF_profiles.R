#Run from the repository root: Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_accessible_SMF_profiles.R
library(data.table)
library(ggplot2)
library(egg)
library(gridExtra)
OUT <- 'output/Fig4_E7_combinedMotifs_insilico/oct4_accessible_profiles'

#Read coverage-weighted profiles from compute_accessible_SMF_profiles.py.
s <- fread(file.path(OUT, 'smoothed_profiles.tsv.gz'))
theme_set(theme_classic() + theme(
  panel.border = element_rect(color = 'black', fill = NA, linewidth = .4),
  axis.line = element_blank(),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 6),
  plot.title = element_text(size = 7),
  plot.subtitle = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank(),
  legend.position = 'none',
  plot.margin = margin(5, 5, 5, 5)))

#Overlay center-only and flanked profiles; shade their protection differences.
make_plot <- function(z, window, title) {
  z <- z[pos >= -window / 2 & pos<window / 2]
  pairs <- dcast(z, spacing + sigma_bp + pos~condition, value.var = 'signal')
  p <- ggplot() +
    geom_ribbon(data = pairs, aes(pos, ymin = A, ymax = pmax(A, BAB)), fill = '#646db3',
      alpha = .22, na.rm = TRUE) +
    geom_ribbon(data = pairs, aes(pos, ymin = pmin(A, BAB), ymax = A), fill = '#454545',
      alpha = .16, na.rm = TRUE) +
    geom_line(data = z[condition == 'BAB'], aes(pos, signal), color = '#646db3', linewidth = .45,
      na.rm = TRUE) +
    geom_line(data = z[condition == 'A'], aes(pos, signal), color = '#454545', linetype = 'dashed',
      linewidth = .4, na.rm = TRUE) +
    scale_x_continuous(limits = c(-window / 2, window / 2), breaks = c(-window / 3, 0, window / 3),
      expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    labs(x = 'Position relative to OCT4 (bp)', y = 'Protection fraction', title = title,
    subtitle = 'Purple: BAB | grey dashed: centre only | shading: difference')
  p
}

#Write the 20-bp and 60-bp examples using the original 4-bp smoothing.
for (cohort_name in 'accessible') for (window in 384) {
  for (partner_name in c('KLF', 'ESRRB')) for (sp in c(20, 60)) {
    z <- s[cohort == cohort_name & partner == partner_name & spacing == sp & sigma_bp == 4]
    p <- make_plot(z, window, paste(partner_name, cohort_name, '|', sp, 'bp spacing | sigma = 4 bp'))
    g <- set_panel_size(p, width = unit(5, 'cm'), height = unit(2.5, 'cm'))
    pdf(file.path(OUT, paste0(partner_name, '_', sp, 'bp_accessible_SMF.pdf')), width = 13 / 2.54,
      height = 7 / 2.54)
    grid.arrange(g, newpage = FALSE)
    dev.off()
  }
}
