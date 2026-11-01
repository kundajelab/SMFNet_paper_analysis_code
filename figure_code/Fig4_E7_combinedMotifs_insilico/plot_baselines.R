#Run from the repository root: Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_baselines.R
library(data.table)
library(ggplot2)
library(egg)
library(gridExtra)
OUT <- "output/Fig4_E7_combinedMotifs_insilico/baselines"

#Read combination medians and matched-background bootstrap statistics.
d <- fread(file.path(OUT, "combination_medians.tsv"))
stats <- fread(file.path(OUT, "baseline_statistics.tsv"))
theme_set(theme_classic() + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
  axis.line = element_blank(),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank(),
  legend.text = element_text(size = 6)))

#Use consistent panel dimensions for the exported PDFs.
save_panel <- function(p, name, width, height) {
  pdf(file.path(OUT, paste0(name, ".pdf")), width = 20 / 2.54, height = 20 / 2.54)
  grid.arrange(set_panel_size(p, width = unit(width, "cm"), height = unit(height, "cm")), newpage = FALSE)
  dev.off()
}

#Scatters use the two spacings displayed in Fig. 4b and E7a,d,f.
for (metric in c("freq", "width", "dispersion")) {
  full <- d[measure == metric]
  z <- full[spacing %in% c(20, 60)]
  z[, highlight := fcase(left == "SOX2", "SOX2", left == "CTCF", "CTCF", default = "Other")]
  if (metric != "freq") z[, highlight := "Other"]
  limits <- if (metric == "freq") c(0, 1) else {
    bounds <- range(c(full$combined_median, full$expected_median))
    bounds + c(-1, 1) * .05 * diff(bounds)
  }
  r <- stats[measure == metric & spacing %in% c(20, 60)]
  p <- ggplot(z, aes(expected_median, combined_median)) + facet_wrap(~spacing, nrow = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "#454545", linewidth = .35) +
    geom_point(data = z[highlight == "Other"], color = "#676767", alpha = .3, size = .75, stroke = 0)
  if (metric == "freq") p <- p +
    geom_point(data = z[highlight != "Other"], aes(color = highlight), alpha = .8, size = 1, stroke = 0) +
    scale_color_manual(values = c(SOX2 = "#667F99", CTCF = "#B08266"), name = "Flank motif")
  p <- p + geom_text(data = r, aes(x = limits[1] + .13 * diff(limits), y = limits[1] + .94 * diff(limits),
    label = paste0("r=", round(pearson_r, 2))), inherit.aes = FALSE, size = 5 / 2.8) +
    coord_cartesian(xlim = limits, ylim = limits) + labs(x = "Baseline", y = "Joint prediction",
      title = metric)
  save_panel(p, paste0(metric, "_baseline_scatter"), 3, 3)
}

#The frequency residual boxplot is shown at 20 bp, grouped by center motif.
z <- d[measure == "freq" & spacing == 20]
z[, residual_pp := 100 * delta]
z[, highlight := fifelse(left %in% c("TEAD", "SOX2", "ETS"), left, "Other")]
levels <- z[, .(med = median(residual_pp)), by = center][order(med), center]
p <- ggplot(z, aes(residual_pp, factor(center, levels = levels))) +
  geom_point(data = z[order(highlight != "Other")], aes(color = highlight),
  position = position_jitter(width = 0, height = .15, seed = 20261021), size = .6, alpha = .8) +
  scale_color_manual(values = c(Other = "#DEDEDE", TEAD = "#D55E00", SOX2 = "#0072B2", ETS = "#009E73"),
  breaks = c("TEAD", "SOX2", "ETS"), name = "Flank motif") +
  geom_boxplot(orientation = "y", fill = NA, outlier.shape = NA, color = "#454545", linewidth = .4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .3) +
  labs(x = "Joint prediction - baseline (pp)", y = "Center motif") + theme(legend.position = "bottom")
save_panel(p, "frequency_residuals_by_center_20bp", 2, 3.5)

#Pearson correlation and signed mean error across all 17 spacings (E7b,e,g).
for (metric in c("freq", "width", "dispersion")) {
  z <- copy(stats[measure == metric])
  for (quantity in c("pearson_r", "mean_error")) {
    if (quantity == "pearson_r") z[, `:=`(estimate = pearson_r, lower = pearson_ci_low,
      upper = pearson_ci_high)]
    else z[, `:=`(estimate = mean_error, lower = mean_error_ci_low, upper = mean_error_ci_high)]
    p <- ggplot(z, aes(factor(spacing), estimate)) + geom_col(fill = "#BCBCBC") +
      geom_linerange(aes(ymin = lower, ymax = upper), linewidth = .25) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
      scale_x_discrete(name = "Motif spacing (bp)", breaks = as.character(seq(20, 100, 20))) +
      labs(y = if (quantity == "pearson_r") "Pearson r" else if (metric == "freq") "Mean error (pp)" else "Mean error (bp)")
    if (quantity == "pearson_r") p <- p + coord_cartesian(ylim = c(0, 1))
    save_panel(p, paste0(metric, "_", quantity), 2, 1.25)
  }
}
