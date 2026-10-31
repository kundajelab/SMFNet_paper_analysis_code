#Run from the repository root: Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_coaccessibility_distance.R
#Run extract_coaccessibility_distances.py first; reads its output table. 
#Writes the distance plot, background-bias table and fitted curves to output/Fig1_E3_Codependency_benchmarking.
library(ggplot2)
library(egg)
library(gridExtra)
library(scales)
library(mgcv)

output_dir <- file.path("output", "Fig1_E3_Codependency_benchmarking")
distance_file <- file.path(output_dir, "coaccessibility_distances.tsv")
if (!file.exists(distance_file)) stop("Missing input: ", distance_file, ". Run extract_coaccessibility_distances.py first.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
#Load only the columns used for fitting, keeping every original pairwise score.
data <- read.delim(distance_file, colClasses = c("NULL", "character", "NULL", "NULL", "numeric", "numeric"))
if (!all(data$label %in% c("ACT", "PRED"))) stop("Unexpected dataset label in ", distance_file)
data$type <- ifelse(data$label == "ACT", "mi_act", "mi_pred")
data$label <- NULL

#Subtract the median score among the longest 20% of distances, separately for each dataset.
bias <- vapply(split(data, data$type), function(x) {
  median(x$mi[x$distance > quantile(x$distance, 0.8)], na.rm = TRUE)
}, numeric(1))
bias_table <- data.frame(type = names(bias), bias = unname(bias))
write.table(bias_table, file.path(output_dir, "coaccessibility_background_bias.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
data$mi_corr <- data$mi - bias[data$type]

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), plot.title = element_text(size = 6), axis.text = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6), strip.text = element_text(size = 6),
  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
))
#Preserve the original plot: ggplot fits the GAM after the pseudo-log scale transformation.
plot <- ggplot(data, aes(distance, mi_corr, color = type)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20), linewidth = 1) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-4),
                    breaks = c(-0.01, -0.001, 0, 0.001, 0.01, 0.1, 1)) +
  scale_x_continuous(limits = c(0, 1024)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#BCBCBC", "#636cb2"))
pdf(file.path(output_dir, "coaccessibility_distance_GAM.pdf"), width = 20 / 2.54, height = 20 / 2.54)
grid.arrange(set_panel_size(plot, height = unit(1.5, "cm"), width = unit(4, "cm")))
dev.off()

#The exported background curves are fitted on the untransformed corrected scores, as before.
fit_table <- do.call(rbind, lapply(split(data, data$type), function(x) {
  fit <- gam(mi_corr ~ s(distance, k = 20), data = x)
  grid <- data.frame(distance = 0:1024)
  prediction <- predict(fit, newdata = grid, se.fit = TRUE)
  data.frame(type = x$type[1], distance = grid$distance, mi_fit = unname(prediction$fit),
             mi_se = unname(prediction$se.fit), lower = unname(prediction$fit - 2 * prediction$se.fit),
             upper = unname(prediction$fit + 2 * prediction$se.fit))
}))
write.table(fit_table, file.path(output_dir, "coaccessibility_background_GAM.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote distance plot, bias and fitted curves to ", output_dir)
