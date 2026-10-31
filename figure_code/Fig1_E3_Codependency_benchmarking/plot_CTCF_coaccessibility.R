#Run from the repository root: Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_CTCF_coaccessibility.R
#Step 3 of 3. Run count_CTCF_motif_states.py and compute_CTCF_coaccessibility.py first.
#Uses input_files/Fig1_E3_Codependency_benchmarking/coaccessibility_background_bias.txt for normalization.
#Region-level four-condition filter (needs reads in all conditions) and 20 bp smoothing radius.
library(ggplot2)
library(dplyr)
library(stringr)
library(egg)
library(gridExtra)
library(scales)

#Configuration

PLOT_DIR <- file.path("output", "Fig1_E3_Codependency_benchmarking", "CTCF")
IN_DIR <- file.path(PLOT_DIR, "matrices")
BIAS_FILE <- file.path("input_files", "Fig1_E3_Codependency_benchmarking", "coaccessibility_background_bias.txt")
if (!file.exists(BIAS_FILE)) stop("Missing normalization input: ", BIAS_FILE)
OUTPUT_FILE <- file.path(PLOT_DIR, "CTCF_aligned_mean_coaccessibility.pdf")

FULL_LEN <- 1024L
TARGET_POS <- 512L
SMOOTH_RADIUS <- 20L

#Require at least this many regions to contribute to a pixel.
#Set to 1 to retain every pixel with at least one contributing region.
MIN_REGIONS_PER_PIXEL <- 1L

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

#Functions

load_sparse_mi_matrix <- function(path) {
  raw <- read.table(path, header = FALSE, sep = "\t", check.names = FALSE)

  col_pos <- as.integer(as.numeric(raw[1, -1]))
  row_pos <- as.integer(as.numeric(raw[-1, 1]))

  mat <- as.matrix(raw[-1, -1])
  mode(mat) <- "numeric"

  list(mat = mat, row_pos = row_pos, col_pos = col_pos)
}

sparse_to_dense <- function(mi_mat, row_pos, col_pos, full_len = 1024L, radius = 20L) {
  mi_full <- matrix(0, nrow = full_len, ncol = full_len)

  mi_count <- matrix(0L, nrow = full_len, ncol = full_len)

  for (i in seq_along(row_pos)) {
    pi <- row_pos[i]

    if (is.na(pi) || pi < 0 || pi >= full_len) {
      next
    }

    x_min <- max(0L, pi - radius)
    x_max <- min(full_len - 1L, pi + radius)
    x_idx <- (x_min:x_max) + 1L

    for (j in seq_along(col_pos)) {
      pj <- col_pos[j]

      if (is.na(pj) || pj < 0 || pj >= full_len) {
        next
      }

      val <- mi_mat[i, j]

      if (is.na(val)) {
        next
      }

      y_min <- max(0L, pj - radius)
      y_max <- min(full_len - 1L, pj + radius)
      y_idx <- (y_min:y_max) + 1L

      mi_full[x_idx, y_idx] <- mi_full[x_idx, y_idx] + val

      mi_count[x_idx, y_idx] <- mi_count[x_idx, y_idx] + 1L
    }
  }

  dense <- mi_full / mi_count
  dense[mi_count == 0] <- NA_real_

  rownames(dense) <- 0:(full_len - 1L)
  colnames(dense) <- 0:(full_len - 1L)

  dense
}

align_matrix_to_motif <- function(mat, motif_pos, full_len = 1024L, target_pos = 512L) {
  if (nrow(mat) != full_len || ncol(mat) != full_len) {
    stop("Matrix dimensions are ", nrow(mat), " x ", ncol(mat), "; expected ", full_len, " x ", full_len)
  }

  aligned <- matrix(NA_real_, nrow = full_len, ncol = full_len)

  #Coordinates are 0-based, while R matrix indices are 1-based.
  original_pos <- 0:(full_len - 1L)
  motif_pos <- as.integer(round(motif_pos))

  #Shift the motif position to TARGET_POS.
  aligned_pos <- original_pos - motif_pos + target_pos

  valid <- (aligned_pos >= 0L & aligned_pos < full_len)

  source_idx <- which(valid)
  destination_idx <- aligned_pos[valid] + 1L

  aligned[destination_idx, destination_idx] <- mat[source_idx, source_idx, drop = FALSE]

  rownames(aligned) <- 0:(full_len - 1L)
  colnames(aligned) <- 0:(full_len - 1L)

  aligned
}

matrix_to_plot_df <- function(mat, type, bound_status, bias, n_regions) {
  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    setNames(c("i", "j", "mi")) %>%
    mutate(i = as.numeric(as.character(i)), j = as.numeric(as.character(j)), mi = mi - bias,
      type = tolower(type), bound_status = bound_status, n_regions_total = n_regions)
}

#Plot theme

theme_set(theme_classic() + theme(legend.text = element_text(size = 6),
      legend.title = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA,
        linewidth = 0.5
      ), panel.background = element_rect(fill = "grey85", color = NA
      ), axis.line = element_blank(),
      axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
      axis.title = element_text(size = 6), plot.title = element_text(size = 7),
      strip.text = element_text(size = 6), strip.background = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "cm")))

#Read bias estimates

cor_vals <- read.delim(BIAS_FILE)

act_bias <- cor_vals %>%
  filter(type == "mi_act") %>%
  pull(bias) %>%
  first()

pred_bias <- cor_vals %>%
  filter(type == "mi_pred") %>%
  pull(bias) %>%
  first()

if (length(act_bias) == 0 || is.na(act_bias)) {
  stop("Could not obtain ACT bias from: ", BIAS_FILE)
}

if (length(pred_bias) == 0 || is.na(pred_bias)) {
  stop("Could not obtain PRED bias from: ", BIAS_FILE)
}

message("ACT bias:  ", act_bias)
message("PRED bias: ", pred_bias)

#Locate and parse matrix files

files <- list.files(IN_DIR, pattern = "__(ACT|PRED)_(bound|unbound)_mi_mat\\.txt$", full.names = TRUE)

if (length(files) == 0) {
  stop("No MI matrix files found in: ", IN_DIR)
}

#Expected example:
#
#region__motif_CTCF__center_512.0__smfpos_510.0
#__ACT_bound_mi_mat.txt
#
filename_pattern <- paste0("^(.*?)", "__motif_(.*?)", "__center_(-?[0-9]+(?:\\.[0-9]+)?)",
  "__smfpos_(-?[0-9]+(?:\\.[0-9]+)?)", "__(ACT|PRED)_(bound|unbound)_mi_mat\\.txt$")

filename_matches <- str_match(basename(files), filename_pattern)

file_meta <- data.frame(file = files, basename = basename(files), region = filename_matches[, 2],
  motif = filename_matches[, 3], motif_center = as.numeric(filename_matches[, 4]),
  motif_smfpos = as.numeric(filename_matches[, 5]), type = filename_matches[, 6],
  bound_status = filename_matches[, 7]) %>%
  filter(!is.na(region), !is.na(motif_smfpos),

    #Exact CTCF only
    motif == "CTCF") %>%
  mutate(condition = paste(type, bound_status, sep = "_")) %>%
  group_by(region) %>%
  filter(all(c("ACT_bound", "ACT_unbound", "PRED_bound", "PRED_unbound"
    ) %in% condition)) %>%
  ungroup() %>%
  arrange(region, type, bound_status)

if (nrow(file_meta) == 0) {
  stop("No exact CTCF regions have all four conditions.")
}

message("CTCF files found: ", nrow(file_meta))

print(file_meta %>%
    count(type, bound_status, name = "n_files"))

#Create sum and count matrices

group_names <- c("ACT_bound", "ACT_unbound", "PRED_bound", "PRED_unbound")

sum_matrices <- setNames(lapply(group_names, function(x) {
      matrix(0, nrow = FULL_LEN, ncol = FULL_LEN)
    }
  ), group_names)

count_matrices <- setNames(lapply(group_names, function(x) {
      matrix(0L, nrow = FULL_LEN, ncol = FULL_LEN)
    }
  ), group_names)

file_counts <- setNames(integer(length(group_names)), group_names)

#Load, densify, align and accumulate every CTCF matrix

for (i in seq_len(nrow(file_meta))) {
  info <- file_meta[i, ]

  message(i, "/", nrow(file_meta), "  Processing: ", info$region, " | ", info$type,
    " | ", info$bound_status, " | motif=", info$motif, " | smfpos=", info$motif_smfpos)

  loaded <- load_sparse_mi_matrix(info$file)

  dense <- sparse_to_dense(mi_mat = loaded$mat, row_pos = loaded$row_pos, col_pos = loaded$col_pos,
    full_len = FULL_LEN, radius = SMOOTH_RADIUS)

  aligned <- align_matrix_to_motif(mat = dense, motif_pos = info$motif_smfpos, full_len = FULL_LEN,
    target_pos = TARGET_POS)

  group <- paste(info$type, info$bound_status, sep = "_")

  valid <- is.finite(aligned)

  if (!any(valid)) {
    warning("No finite values after alignment for: ", info$file)

    next
  }

  sum_matrices[[group]][valid] <- sum_matrices[[group]][valid] + aligned[valid]

  count_matrices[[group]][valid] <- count_matrices[[group]][valid] + 1L

  file_counts[[group]] <- file_counts[[group]] + 1L
}

message("Successfully processed files:")

print(file_counts)

#Calculate mean matrices across all aligned CTCF instances

mean_matrices <- setNames(vector("list", length(group_names)), group_names)

for (group in group_names) {
  if (file_counts[[group]] == 0) {
    warning("No usable files for group: ", group)

    mean_matrices[[group]] <- matrix(NA_real_, nrow = FULL_LEN, ncol = FULL_LEN)

    next
  }

  mean_mat <- sum_matrices[[group]] /
    count_matrices[[group]]

  #Leave pixels as NA when too few regions contributed.
  mean_mat[
    count_matrices[[group]] < MIN_REGIONS_PER_PIXEL
  ] <- NA_real_

  rownames(mean_mat) <- 0:(FULL_LEN - 1L)
  colnames(mean_mat) <- 0:(FULL_LEN - 1L)

  mean_matrices[[group]] <- mean_mat
}

#Create combined plotting table

plot_list <- list()

for (type in c("ACT", "PRED")) {
  for (bound_status in c("bound", "unbound")) {
    group <- paste(type, bound_status, sep = "_")

    mat <- mean_matrices[[group]]

    bias <- if (type == "ACT") {
      act_bias
    } else {
      pred_bias
    }

    plot_list[[group]] <- matrix_to_plot_df(mat = mat, type = type, bound_status = bound_status,
      bias = bias, n_regions = file_counts[[group]])
  }
}

df_plot <- bind_rows(plot_list) %>%
  mutate(type = factor(type, levels = c("act", "pred")
    ), bound_status = factor(bound_status, levels = c("bound", "unbound")
    ))

if (!any(is.finite(df_plot$mi))) {
  stop("No finite values remained after alignment and averaging.")
}

breaks <- c(-1, -0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.1, 1)

#Combined CTCF-aligned plot

ac_heat <- df_plot %>%
  ggplot(aes(x = as.numeric(i), y = as.numeric(j), fill = mi)) + geom_raster() + facet_grid(type ~ bound_status,
    labeller = labeller(type = c(act = "Experiment", pred = "Model"),
                        bound_status = c(bound = "Bound", unbound = "Unbound"))) +
  geom_vline(xintercept = TARGET_POS, color = "black", linewidth = 0.15,
    linetype = "dashed") + geom_hline(yintercept = TARGET_POS,
    color = "black", linewidth = 0.15, linetype = "dashed") + scale_x_continuous(breaks = c(0, 256,
      TARGET_POS, 768, 1023
    ), labels = c("-512", "-256", "CTCF", "+256", "+511"
    ), limits = c(0, 1023), expand = c(0, 0)) + scale_y_continuous(breaks = c(0, 256,
      TARGET_POS, 768, 1023
    ), labels = c("-512", "-256", "CTCF", "+256", "+511"
    ), limits = c(0, 1023), expand = c(0, 0)) +
  scale_fill_gradient2(low = "#008837", mid = "white", high = "#7b3294",
    midpoint = 0, na.value = "grey85", trans = pseudo_log_trans(sigma = 1e-4
    ), limits = c(-1, 1
    ), breaks = breaks, oob = squish) +
  coord_fixed() + labs(title = paste0("Mean CTCF-aligned co-accessibility",
      " (radius = ", SMOOTH_RADIUS, ")"
    ), x = "Position relative to CTCF motif", y = "Position relative to CTCF motif", fill = "Mean score"
  )

pdf(OUTPUT_FILE, width = 20 / 2.54, height = 20 / 2.54)

fixed_plot <- set_panel_size(ac_heat, height = unit(2, "cm"), width = unit(2, "cm"))

grid.arrange(fixed_plot)

dev.off()

message("Done")
message("Wrote combined plot: ", OUTPUT_FILE)
