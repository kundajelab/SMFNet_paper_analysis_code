library(ggplot2)
library(egg)
library(gridExtra)
library(data.table)
library(viridis)
library(ggrastr)

#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_SMF-ONT_example.R
input_dir <- file.path("input_files", "Fig1_E1_ONT_benchmarking")
output_dir <- file.path("output", "Fig1_E1_ONT_benchmarking")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(),
  axis.text = element_text(size = 6),
  title = element_text(size = 6),
  axis.title = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank()
))

save_panel <- function(plot, filename, height, width) {
  pdf(file.path(output_dir, filename), width = 20 / 2.54, height = 20 / 2.54)
  on.exit(dev.off())
  panel <- set_panel_size(plot, height = grid::unit(height, "cm"),
                          width = grid::unit(width, "cm"))
  grid.arrange(panel)
}

length_files <- file.path(input_dir, paste0("SMF-ONT_R", 1:4, "_read_lengths.txt"))
stopifnot(all(file.exists(length_files)))
read_lengths <- rbindlist(lapply(length_files, function(path) {
  data <- fread(path, header = FALSE, col.names = "read_length")
  data[, sample := sub("_read_lengths\\.txt$", "", basename(path))]
  data
}))

#N50 is the read length at which half the sequenced bases are accounted for.
calc_n50 <- function(x) {
  x <- sort(as.numeric(x), decreasing = TRUE)
  x[which(cumsum(x) >= sum(x) / 2)[1]]
}
genome_size <- 2.7e9
n50_df <- read_lengths[, .(
  N50 = calc_n50(read_length),
  read_n = .N,
  total_bases = sum(read_length),
  coverage = sum(read_length) / genome_size,
  median = median(read_length),
  ymax = max(read_length)
), by = sample]

length_plot <- ggplot(read_lengths, aes(x = sample, y = read_length)) +
  geom_violin(fill = "#DEDEDE") +
  stat_summary(fun = calc_n50, geom = "crossbar", linewidth = 0.5, width = 0.25) +
  geom_text(data = n50_df,
            aes(x = sample, y = 40000, label = paste0("N50 = ", round(N50))),
            inherit.aes = FALSE, size = 6 / 2.8) +
  scale_y_continuous(limits = c(0, 50000))
save_panel(length_plot, "length_qc_ont.pdf", 2.5, 3)

#Use sites with complete data across all replicates and the short-read reference.
cor_df <- na.omit(fread(file.path(
  input_dir, "SMF-ONT_TKO_avg_meth_vs_all_minCov_20.tsv.gz"
)))
stopifnot(ncol(cor_df) == 11L)
setnames(cor_df, c("site_index", "chr", "start", "end", "strand", "context",
                   "R1", "R2", "R3", "R4", "Sonmezer"))

correlation_plot <- function(data, x, y, bins, fill_limits = NULL) {
  correlation <- cor(data[[x]], data[[y]])
  ggplot(data, aes(x = 1 - .data[[x]], y = 1 - .data[[y]])) +
    rasterize(geom_hex(bins = bins), dpi = 1200) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_viridis(option = "A", trans = "log10", limits = fill_limits,
                       oob = scales::squish) +
    annotate("text", label = paste0("r=", round(correlation, 2)),
             size = 5 / 2.8, x = 0.15, y = 0.9)
}

#Keep both orientations of each replicate pair, as in the original plots.
replicates <- paste0("R", 1:4)
for (rep_a in replicates) {
  for (rep_b in setdiff(replicates, rep_a)) {
    plot <- correlation_plot(cor_df, rep_a, rep_b, bins = 200,
                             fill_limits = c(1, 2000))
    save_panel(plot, paste0("cor_reps_", rep_a, "_", rep_b, ".pdf"), 1.5, 1.5)
  }
}

#Compare the equally weighted ONT replicate mean with the Sonmezer reference.
cor_df[, ont := (R1 + R2 + R3 + R4) / 4]
reference_plot <- correlation_plot(cor_df, "ont", "Sonmezer", bins = 500)
save_panel(reference_plot, "cor_sonmezer_ont.pdf", 2, 2)
