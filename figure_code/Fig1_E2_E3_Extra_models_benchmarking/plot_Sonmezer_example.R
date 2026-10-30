#Run from the repository root: Rscript figure_code/Fig1_E2_E3_Extra_models_benchmarking/plot_Sonmezer_example.R
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(egg)
library(gridExtra)
library(hdf5r)

input_dir <- file.path("input_files", "Fig1_E2_E3_Extra_models_benchmarking")
output_dir <- file.path("output", "Fig1_E2_E3_Extra_models_benchmarking")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
#Read smfnet.samples v2 matrices as molecules by positions.
#hdf5r reverses the stored dimensions; positions retain their original convention.
sample_index <- function(h5) {
  stopifnot(h5attr(h5, "schema") == "smfnet.samples")
  #Map genomic region IDs to the numbered groups stored in the H5 file.
  ids <- as.character(h5[["region_ids"]][])
  setNames(sprintf("regions/%08d", seq_along(ids) - 1L), ids)
}
read_sample_region <- function(h5, index, region, role) {
  stopifnot(region %in% names(index))
  group <- h5[[index[[region]]]]
  stopifnot(h5attr(group, "region_id") == region)
  #The role selects observed or predicted reads within the region group.
  mat <- t(group[[paste0(role, "/smf")]]$read(drop = FALSE))
  pos <- as.numeric(group[[paste0(role, "/positions")]][])
  stopifnot(ncol(mat) == length(pos))
  list(mat = mat, pos = pos)
}

#Reuse the motif annotations shared with the ONT benchmarking figures.
tfbs_file <- file.path("input_files", "Fig1_E1_ONT_benchmarking", "SMFNet_annotated_motif_hits.tsv")

tfbs_df <- read.delim(tfbs_file) %>%
  select(chr, start, end, pattern, motif)

theme_set(
  theme_classic() +
    theme(
      legend.text = element_text(size = 6),
      panel.border = element_rect(
        color = "black", fill = NA, linewidth = 0.5
      ),
      axis.line = element_blank(),
      axis.text.x = element_text(size = 6),
      plot.title = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      title = element_text(size = 6),
      legend.position = "none",
      axis.title = element_text(size = 6),
      strip.text = element_text(size = 6),
      strip.background = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
)

#The original example uses a 1-based observed region label. The publication
#files index both datasets by the corresponding 0-based prediction region ID.
region <- "chr1:43137426-43139473"
sample_region <- "chr1:43137425-43139473"

region_parts <- strsplit(region, "[:-]")[[1]]
zoom_region <- paste0(region_parts[1], ":", as.numeric(region_parts[2]) + 512, "-", as.numeric(region_parts[3]) - 512)
zoom_region_safe <- gsub("[:-]", "_", zoom_region)

#Express overlapping motif hits relative to the example region.
get_tfbs_for_region <- function(region) {
  coordinates <- strsplit(region, "[:-]")[[1]]
  chromosome <- coordinates[1]
  region_start <- as.numeric(coordinates[2])
  region_end <- as.numeric(coordinates[3])
  #Preserve the original inclusive overlap and legacy site coordinates.
  tfbs_df %>%
    filter(chr == chromosome, start <= region_end, end >= region_start) %>%
    mutate(rel_start = start - region_start,
           rel_end = end - region_start,
           rel_center = (rel_start + rel_end) / 2) %>%
    filter(rel_end >= 0, rel_start < 2048) %>%
    mutate(rel_start = pmax(rel_start, 0), rel_end = pmin(rel_end, 2048))
}

pred_handle <- H5File$new(file.path(input_dir, "mESC_Sonmezer_testRegions_SMFNet_reads.h5"), mode = "r")
act_handle <- H5File$new(file.path(input_dir, "mESC_Sonmezer_testRegions_Experiment_reads.h5"), mode = "r")
pred_index <- sample_index(pred_handle)
act_index <- sample_index(act_handle)

#Load all reads for the selected example; shift observed sites to the prediction window.
pred <- read_sample_region(pred_handle, pred_index, sample_region, "predicted")
act <- read_sample_region(act_handle, act_index, sample_region, "observed")
act$pos <- act$pos - 512

tfbs_region <- get_tfbs_for_region(region) %>%
  mutate(rel_start = rel_start - 512, rel_end = rel_end - 512, rel_center = rel_center - 512) %>%
  filter(rel_end >= 0 & rel_end <= 1024)

#Convert matrices to one row per molecule and site for plotting.
pred_long <- as.data.frame(pred$mat)
colnames(pred_long) <- pred$pos

pred_long <- pred_long %>%
  mutate(read_id = row_number()) %>%
  pivot_longer(
    -read_id,
    names_to = "pos",
    values_to = "value"
  ) %>%
  mutate(
    pos = as.numeric(pos),
    type = "pred"
  )

act_long <- as.data.frame(act$mat)
colnames(act_long) <- act$pos

act_long <- act_long %>%
  mutate(read_id = row_number()) %>%
  pivot_longer(
    -read_id,
    names_to = "pos",
    values_to = "value"
  ) %>%
  mutate(
    pos = as.numeric(pos),
    type = "act"
  )

sm_long <- bind_rows(pred_long, act_long) %>%
  filter(pos >= 0, pos < 1024) %>%
  mutate(molecule_id = paste(type, read_id, sep = "_"))

#Cluster observed and predicted molecules together to obtain a shared ordering.
cluster_mat <- sm_long %>%
  select(molecule_id, pos, value) %>%
  pivot_wider(
    names_from = pos,
    values_from = value
  ) %>%
  column_to_rownames("molecule_id") %>%
  as.matrix()

#Preserve the original clustering: only NA entries are excluded from overlap.
pairwise_smf_dist <- function(mat, min_overlap = 10) {

  mat <- as.matrix(mat)

  #Observed positions
  obs <- !is.na(mat)

  #Binary states, missing -> 0 only for matrix operations
  x <- mat
  x[!obs] <- 0

  obs <- obs * 1

  #Number of jointly observed positions
  overlap <- tcrossprod(obs)

  #Number of shared 1/1 positions
  both_one <- tcrossprod(x)

  #Number of 1s in read i at positions observed in read j
  one_i <- x %*% t(obs)

  #Binary disagreement:
  #(1,0) + (0,1)
  disagreement <- one_i + t(one_i) - 2 * both_one

  #Fraction disagreement over jointly observed positions
  d <- disagreement / overlap

  #Insufficient overlap -> maximal distance
  d[overlap < min_overlap] <- 1
  d[!is.finite(d)] <- 1

  diag(d) <- 0

  rownames(d) <- rownames(mat)
  colnames(d) <- rownames(mat)

  d
}

dmat <- pairwise_smf_dist(
  cluster_mat,
  min_overlap = 10
)

dmat[is.na(dmat)] <- 1

hc <- hclust(
  as.dist(dmat),
  method = "average"
)

cluster_order <- tibble(
  molecule_id = rownames(cluster_mat)[hc$order],
  global_order = seq_along(hc$order)
)

sm_long_plot <- sm_long %>%
  left_join(
    cluster_order,
    by = "molecule_id"
  ) %>%
  group_by(type) %>%
  mutate(
    order = dense_rank(global_order)
  ) %>%
  ungroup()

#Widen each site marker for visibility in the single-molecule panels.
sm_long_plot_raster <- sm_long_plot %>%
  tidyr::crossing(dx = -2.5:2.5) %>%
  mutate(pos = pos + dx)

sm_plot <- sm_long_plot_raster %>%
  ggplot(aes(x = pos, y = order, fill = as.factor(value))) +
  facet_wrap(~type, scales = "free") +
  geom_raster() +
  scale_fill_manual(
    values = c("1" = "#DEDEDE",
               "0" = "#000000",
               "-1" = "#FFFFFF")
  ) +
  coord_cartesian(xlim = c(0, 1024), expand = FALSE) +
  scale_x_continuous(
    breaks = c(0, 256, 512, 768, 1024)
  ) +
  scale_y_continuous(
    name = "Single molecules",
    expand = c(0, 0)
  ) +
  labs(
    x = "Relative position (bp)",
    title = paste0(zoom_region)
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

pdf(file.path(output_dir, paste0("Sonmezer_example_", zoom_region_safe, "_sm.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(sm_plot, height = unit(3, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
dev.off()

#Compute site averages from all reads, excluding missing calls (-1).
get_average_from_matrix <- function(mat, pos, type) {

  mat[mat == -1] <- NA

  tibble(
    smf_pos = pos,
    avg_signal = colMeans(mat, na.rm = TRUE),
    type = type
  )
}

act_region_avg <- get_average_from_matrix(act$mat, act$pos, "act")
act_region_pred <- get_average_from_matrix(pred$mat, pred$pos, "pred")

avg_combi <- bind_rows(act_region_avg, act_region_pred)

p_avg <- avg_combi %>%
  ggplot(aes(x = smf_pos, y = avg_signal, color = type)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 0.4) +
  geom_rect(
    data = tfbs_region,
    inherit.aes = FALSE,
    aes(
      xmin = rel_start,
      xmax = rel_end,
      ymin = 0,
      ymax = 0.08
    ),
    fill = "#BCBCBC"
  ) +

  geom_text(
    data = tfbs_region,
    inherit.aes = FALSE,
    aes(
      x = rel_center,
      y = 0.12,
      label = motif
    ),
    size = 2
  ) +
  scale_color_manual(
    values = c("act" = "#bcbcbc",
      "pred" = "#636ab1")) +
  scale_y_continuous(limits = c(0, 1), name = "Mean SMF signal") +
  scale_x_continuous(limits = c(0, 1024), breaks = c(0, 256, 512, 768, 1024), name = "Position") +
  theme(legend.position = "right") +
  ggtitle(zoom_region)

pdf(file.path(output_dir, paste0("Sonmezer_example_", zoom_region_safe, "_avg.pdf")),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(p_avg, height = unit(3, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
dev.off()

#Close files
pred_handle$close_all()
act_handle$close_all()
