#Run from the repository root: Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_coaccessibility_matrices.R
#Run compute_coaccessibility.py and smooth_coaccessibility_matrices.py first.
#Reads generated smoothed matrices and the supplied background-bias TXT file.
library(ggplot2)
library(egg)
library(gridExtra)
library(scales)

pipeline_dir <- file.path("output", "Fig1_E3_Codependency_benchmarking")
output_dir <- file.path(pipeline_dir, "coaccessibility_matrices_plots")
matrix_dir <- file.path(pipeline_dir, "coaccessibility_smoothed_matrices")
bias_file <- file.path("input_files", "Fig1_E3_Codependency_benchmarking", "coaccessibility_background_bias.txt")
if (!file.exists(bias_file)) stop("Missing normalization input: ", bias_file)
bias_table <- read.delim(bias_file)
bias <- c(ACT = bias_table$bias[match("mi_act", bias_table$type)],
          PRED = bias_table$bias[match("mi_pred", bias_table$type)])
if (any(!is.finite(bias))) stop("Missing experiment/model bias in ", bias_file)
files <- list.files(matrix_dir, pattern = "_(ACT|PRED)_mi_1024\\.txt$", full.names = TRUE)
if (!length(files)) stop("No smoothed matrices in ", matrix_dir, ". Run smooth_coaccessibility_matrices.py first.")
regions <- unique(sub("_(ACT|PRED)_mi_1024\\.txt$", "", basename(files)))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), axis.text = element_text(size = 6), plot.title = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6), strip.text = element_text(size = 6),
  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
))

read_plot_matrix <- function(path, label) {
  raw <- read.table(path, header = FALSE, sep = "\t", check.names = FALSE)
  matrix <- as.matrix(raw[-1, -1])
  mode(matrix) <- "numeric"
  dimnames(matrix) <- list(as.numeric(raw[-1, 1]), as.numeric(raw[1, -1]))
  data <- as.data.frame(as.table(matrix))
  names(data) <- c("i", "j", "mi")
  data$mi <- data$mi - bias[[label]]
  data$type <- tolower(label)
  data
}

#Read and plot one experiment/model pair at a time, retaining the original axes and colour scale.
for (region in regions) {
  paths <- file.path(matrix_dir, paste0(region, c("_ACT_mi_1024.txt", "_PRED_mi_1024.txt")))
  if (!all(file.exists(paths))) stop("Missing experiment/model matrix for ", region)
  data <- rbind(read_plot_matrix(paths[1], "ACT"), read_plot_matrix(paths[2], "PRED"))
  if (!any(is.finite(data$mi))) stop("No finite scores for ", region)
  limit <- max(abs(data$mi), na.rm = TRUE)
  if (limit == 0) limit <- 1
  plot <- ggplot(data, aes(as.numeric(i), as.numeric(j), fill = mi)) +
    geom_raster() + facet_wrap(~type) +
    scale_y_continuous(breaks = c(0, 500, 1000)) +
    scale_x_continuous(breaks = c(0, 500, 1000)) +
    scale_fill_gradient2(low = "#008837", mid = "white", high = "#7b3294",
      trans = pseudo_log_trans(sigma = 1e-4), limits = c(-limit, limit),
      breaks = c(-1, -0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.1, 1)) +
    coord_fixed()
  path <- file.path(output_dir, paste0(region, "_coaccessibility.pdf"))
  pdf(path, width = 20 / 2.54, height = 20 / 2.54)
  grid.arrange(set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm")))
  dev.off()
  message("Wrote: ", path)
}
