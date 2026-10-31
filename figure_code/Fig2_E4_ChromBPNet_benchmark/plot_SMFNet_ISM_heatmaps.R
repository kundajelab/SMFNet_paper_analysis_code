#Run from the repository root: Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_SMFNet_ISM_heatmaps.R
#Plot every raw ISM matrix over its full output window, without score-based selection.
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(egg)
library(gridExtra)
library(hdf5r)
library(GenomicRanges)
library(IRanges)

input_dir <- file.path("input_files", "Fig2_E4_ChromBPNet_benchmark")
output_dir <- file.path("output", "Fig2_E4_ChromBPNet_benchmark", "SMFNet_ISM_heatmaps")
ism_file <- file.path(input_dir, "SMFNet_endogenous_ISM.h5")
reads_file <- file.path(input_dir, "mESC_SMF-ONT_testRegions_SMFNet_reads.h5")
tfbs_file <- file.path(input_dir,
                       "SMFNet_annotated_motif_hits_for_codependency_analysis.tsv")

theme_set(theme_classic() + theme(
  text = element_text(size = 6), axis.text = element_text(size = 6),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), strip.background = element_blank(),
  plot.margin = unit(c(0, 0, 0, 0), "cm")
))

save_panel <- function(plot, filename, height, width = 4) {
  pdf(file.path(output_dir, filename), width = 20 / 2.54, height = 20 / 2.54)
  on.exit(dev.off())
  grid.arrange(set_panel_size(plot, height = unit(height, "cm"), width = unit(width, "cm")))
}

#Read the packaged samples as molecules by positions, preserving the stored calls.
read_sample_region <- function(h5, index, region) {
  group <- h5[[index[[region]]]]
  mat <- t(group[["predicted/smf"]]$read(drop = FALSE))
  pos <- as.numeric(group[["predicted/positions"]][])
  stopifnot(ncol(mat) == length(pos))
  list(mat = mat, pos = pos)
}

average_region <- function(data, region_id, type) {
  data %>% filter(region == region_id) %>%
    mutate(smf_pos = str_split(smf_pos, ","), avg_signal = str_split(avg_signal, ",")) %>%
    unnest(c(smf_pos, avg_signal)) %>%
    mutate(smf_pos = as.numeric(smf_pos), avg_signal = as.numeric(avg_signal), type = type)
}

main <- function() {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tfbs <- read.delim(tfbs_file)
  motifs <- GRanges(tfbs$chr, IRanges(tfbs$start, tfbs$end), motif = tfbs$motif)
  pred_avg <- read.delim(file.path(input_dir, "avg_TestRegions_SMFNet.tsv"))
  act_avg <- read.delim(file.path(input_dir, "avg_TestRegions_Experiment.tsv"),
                       header = FALSE, col.names = c("region", "smf_pos", "avg_signal"))
  ism_h5 <- H5File$new(ism_file, mode = "r")
  on.exit(ism_h5$close_all(), add = TRUE)
  reads_h5 <- H5File$new(reads_file, mode = "r")
  on.exit(reads_h5$close_all(), add = TRUE)
  stopifnot(h5attr(reads_h5, "schema") == "smfnet.samples")
  read_ids <- as.character(reads_h5[["region_ids"]][])
  read_index <- setNames(sprintf("regions/%08d", seq_along(read_ids) - 1L), read_ids)
  regions <- names(ism_h5[["ism_mat"]])

  for (i in seq_along(regions)) {
    rid <- regions[i]
    message("[", i, "/", length(regions), "] ", rid)
    coords <- str_match(rid, "^(chr[^_]+)_([0-9]+)-([0-9]+)$")
    stopifnot(!anyNA(coords))
    #The 1024-bp model output starts 512 bp into the original 2048-bp region.
    start <- as.numeric(coords[1, 3]) + 512
    ism <- ism_h5[["ism_mat"]][[rid]]$read(drop = FALSE)
    stopifnot(nrow(ism) == 1024L, ncol(ism) == 1024L)
    window <- GRanges(coords[1, 2], IRanges(start, start + 1023))
    hits <- findOverlaps(motifs, window)
    annotation <- as.data.frame(motifs[queryHits(hits)]) %>%
      mutate(rel_start = start - .env$start, rel_end = end - .env$start) %>%
      filter(rel_start >= 0, rel_end < 1024)

    #Use raw ism_mat values directly; neither axes nor colors use scaled scores.
    heat <- expand.grid(row = 0:1023, col = 0:1023)
    heat$value <- as.vector(ism)
    limit <- max(abs(heat$value), na.rm = TRUE) * 1.1
    if (!is.finite(limit) || limit == 0) limit <- 1
    heatmap <- ggplot(heat, aes(row, col, fill = value)) +
      geom_raster() +
      scale_fill_gradient2(low = "#008837", mid = "white", high = "#7b3294",
                           midpoint = 0, limits = c(-limit, limit)) +
      coord_equal() +
      labs(x = "Mutated position (bp)", y = "Output position (bp)", fill = "Raw ISM") +
      geom_rect(data = annotation, inherit.aes = FALSE,
                aes(xmin = rel_start, xmax = rel_end, ymin = 5, ymax = 20), fill = "#BCBCBC") +
      geom_text(data = annotation, inherit.aes = FALSE,
                aes(x = (rel_start + rel_end) / 2, y = 25, label = motif), size = 4 / 2.8) +
      geom_rect(data = annotation, inherit.aes = FALSE,
                aes(ymin = rel_start, ymax = rel_end, xmin = 5, xmax = 20), fill = "#BCBCBC") +
      geom_text(data = annotation, inherit.aes = FALSE,
                aes(y = (rel_start + rel_end) / 2, x = 25, label = motif), size = 4 / 2.8)
    save_panel(heatmap, paste0(rid, "_heat.pdf"), height = 4)

    #Companion profiles use the same full window and the existing sampled reads.
    region_id <- sub("_", ":", rid, fixed = TRUE)
    if (!region_id %in% names(read_index)) {
      message("No sampled reads for ", region_id, "; heatmap was written.")
      next
    }
    sample <- read_sample_region(reads_h5, read_index, region_id)
    keep <- sample$pos >= 0 & sample$pos < 1024
    mat <- sample$mat[, keep, drop = FALSE]
    pos <- sample$pos[keep]
    if (!length(pos) || !nrow(mat)) next
    averages <- bind_rows(average_region(act_avg, region_id, "act"),
                          average_region(pred_avg, region_id, "pred")) %>%
      filter(smf_pos %in% pos)
    avg_plot <- ggplot(averages, aes(smf_pos, avg_signal, color = type)) +
      geom_line() + geom_point(size = 0.5) +
      geom_rect(data = annotation, inherit.aes = FALSE,
                aes(xmin = rel_start, xmax = rel_end, ymin = 0, ymax = 0.1), fill = "#BCBCBC") +
      geom_text(data = annotation, inherit.aes = FALSE,
                aes(x = (rel_start + rel_end) / 2, y = 0.15, label = motif), size = 4 / 2.8) +
      scale_color_manual(values = c(act = "#CDCDCD", pred = "#636ab1")) +
      scale_y_continuous(limits = c(0, 1), name = "Mean SMF signal") +
      scale_x_continuous(limits = c(0, 1024), name = "Relative position (bp)")
    save_panel(avg_plot, paste0(rid, "_avg.pdf"), height = 2)

    #Retain the original Manhattan-distance, average-linkage read ordering.
    order <- if (nrow(mat) > 1) hclust(dist(mat, method = "manhattan"), method = "average")$order else 1L
    ordered_reads <- as.data.frame(mat)
    names(ordered_reads) <- as.character(pos)
    ordered_reads$read_id <- paste0("read_id_", seq_len(nrow(mat)))
    ordered_reads$order <- match(seq_len(nrow(mat)), order)
    sm_plot <- ordered_reads %>%
      pivot_longer(-c(read_id, order), values_to = "state", names_to = "smf_pos") %>%
      mutate(smf_pos = as.numeric(smf_pos)) %>%
      ggplot(aes(smf_pos, order, fill = factor(state))) + geom_tile() +
      scale_fill_manual(values = c("-1" = "white", "0" = "#DEDEDE", "1" = "#000000")) +
      coord_cartesian(xlim = c(0, 1024)) +
      labs(x = "Relative position (bp)", y = "Single molecules", fill = "State") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    save_panel(sm_plot, paste0(rid, "_sm.pdf"), height = 2)
    write.table(ordered_reads, file.path(output_dir, paste0(rid, "_reads.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

main()
