#Run from the repository root: Rscript figure_code/Fig1_E2_E3_Extra_models_benchmarking/plot_fiber_regions.R
#Optional arguments select H5 region group IDs, for example: 00000000 00000001.
#Without arguments, plot the original reproducible sample of 200 regions.
library(ggplot2)
library(dplyr)
library(tibble)
library(hdf5r)
library(ggrastr)
library(egg)
library(gridExtra)

input_dir <- file.path("input_files", "Fig1_E2_E3_Extra_models_benchmarking")
output_dir <- file.path("output", "Fig1_E2_E3_Extra_models_benchmarking", "fiber_example_regions")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

CONDITIONS <- c("observed", "predicted")
READ_COLORS <- c("0" = "#2c7fb8", "1" = "#ffffcc")

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), axis.text = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6),
  strip.text = element_text(size = 6), strip.background = element_blank()
))

#Read the existing samples-v2 layout; hdf5r reverses matrix dimensions.
read_condition <- function(group) {
  pos <- as.numeric(group[["positions"]][])
  mat <- t(group[["smf"]]$read(drop = FALSE))
  stopifnot(length(pos) > 0, ncol(mat) == length(pos), nrow(mat) > 0,
            !anyDuplicated(pos))
  order <- order(pos)
  list(pos = pos[order], mat = mat[, order, drop = FALSE])
}

#Align sites before pooling reads; unavailable sites and -1 calls remain missing.
pool_conditions <- function(data) {
  positions <- sort(unique(unlist(lapply(data, function(x) x$pos))))
  matrices <- lapply(data, function(x) {
    mat <- matrix(NA_real_, nrow(x$mat), length(positions))
    mat[, match(x$pos, positions)] <- x$mat
    mat[mat == -1] <- NA_real_
    mat
  })
  list(mat = do.call(rbind, matrices), pos = positions,
       condition = rep(names(data), vapply(matrices, nrow, integer(1))))
}

#Cluster 32 bp bin averages jointly, then place sparsely covered reads last.
order_molecules <- function(pooled, positions) {
  bins <- split(seq_along(positions), floor((positions - min(positions)) / 32))
  binned <- do.call(cbind, lapply(bins, function(cols) {
    mat <- pooled[, cols, drop = FALSE]
    signal <- rowMeans(mat, na.rm = TRUE)
    signal[rowSums(!is.na(mat)) < 2] <- NA_real_
    signal
  }))
  coverage <- rowMeans(is.finite(binned))
  zero_fraction <- rowMeans(pooled == 0, na.rm = TRUE)
  zero_fraction[!is.finite(zero_fraction)] <- 0
  good <- which(coverage >= 0.2)
  sparse <- which(coverage < 0.2)
  ordered_good <- good

  #A dendrogram requires at least two sufficiently covered molecules.
  if (length(good) >= 2) {
    mat <- binned[good, , drop = FALSE]
    distances <- matrix(0, nrow(mat), nrow(mat))
    for (i in seq_len(nrow(mat) - 1L)) {
      for (j in (i + 1L):nrow(mat)) {
        shared <- is.finite(mat[i, ]) & is.finite(mat[j, ])
        distances[i, j] <- distances[j, i] <- if (sum(shared) >= 3) {
          mean(abs(mat[i, shared] - mat[j, shared]))
        } else 1
      }
    }
    dend <- as.dendrogram(hclust(as.dist(distances), method = "ward.D2"))
    dend <- reorder(dend, wts = -zero_fraction[good], agglo.FUN = mean)
    ordered_good <- good[order.dendrogram(dend)]
  }
  c(ordered_good, sparse[order(coverage[sparse], decreasing = TRUE)])
}

save_panel <- function(plot, filename, height, width) {
  pdf(file.path(output_dir, filename), width = 20 / 2.54, height = 20 / 2.54)
  on.exit(dev.off())
  panel <- set_panel_size(plot, height = grid::unit(height, "cm"), width = grid::unit(width, "cm"))
  grid.arrange(panel)
}

main <- function(selected_ids = commandArgs(trailingOnly = TRUE)) {
  h5 <- H5File$new(file.path(input_dir, "GM12878_Fiber-seq_testRegions_samples.h5"), mode = "r")
  on.exit(h5$close_all())
  stopifnot(h5attr(h5, "schema") == "smfnet.samples")
  regions <- h5[["regions"]]
  region_labels <- as.character(h5[["region_ids"]][])
  ids <- sprintf("%08d", seq_along(region_labels) - 1L)
  if (length(selected_ids)) {
    stopifnot(all(selected_ids %in% ids))
    ids <- unique(selected_ids)
  } else {
    set.seed(123)
    ids <- sample(ids, min(200, length(ids)))
  }

  for (id in ids) {
    region_label <- region_labels[as.integer(id) + 1L]
    message("Processing ", id, ": ", region_label)
    data <- setNames(lapply(CONDITIONS, function(cond) {
      read_condition(regions[[id]][[cond]])
    }), CONDITIONS)

    #Average all available calls at each site, excluding missing values.
    avg_df <- bind_rows(lapply(CONDITIONS, function(cond) {
      sample <- data[[cond]]
      sample$mat[sample$mat == -1] <- NA_real_
      tibble(condition = cond, pos = sample$pos, signal = colMeans(sample$mat, na.rm = TRUE))
    }))
    avg_plot <- ggplot(avg_df, aes(pos, signal, color = condition)) +
      geom_line(linewidth = 0.4) +
      scale_color_manual(values = c(observed = "#BCBCBC", predicted = "#646db3")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(breaks = c(0, 256, 512, 768, 1024)) +
      labs(x = "Position (bp)", y = "Mean Fiber-seq signal", color = NULL, title = region_label) +
      theme(legend.position = "bottom")

    pooled <- pool_conditions(data)
    leaf_order <- order_molecules(pooled$mat, pooled$pos)
    #Use the shared dendrogram order within each condition's panel.
    sm_df <- bind_rows(lapply(CONDITIONS, function(cond) {
      idx <- leaf_order[pooled$condition[leaf_order] == cond]
      columns <- match(data[[cond]]$pos, pooled$pos)
      mat <- pooled$mat[idx, columns, drop = FALSE]
      n <- nrow(mat)
      n_pos <- ncol(mat)
      tibble(condition = cond, pos = rep(data[[cond]]$pos, times = n),
             y = rep(1 - (seq_len(n) - 0.5) / n, each = n_pos), height = 1 / n,
             state = factor(as.vector(t(mat)), levels = c(0, 1)))
    })) %>% mutate(condition = factor(condition, levels = CONDITIONS))

    sm_plot <- ggplot(sm_df, aes(pos, y, fill = state)) +
      rasterize(geom_tile(aes(height = height), width = 5), dpi = 1200) +
      facet_wrap(~condition, nrow = 1) +
      scale_fill_manual(values = READ_COLORS, na.value = "white", drop = FALSE) +
      scale_x_continuous(expand = c(0, 0), breaks = c(0, 256, 512, 768, 1024)) +
      #Clip at the panel boundary without dropping tiles due to rounding.
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Position (bp)", y = "Single molecules", title = region_label) +
      theme(legend.position = "none")

    save_panel(avg_plot, paste0("Fiber-seq_", id, "_avg.pdf"), 3, 6)
    save_panel(sm_plot, paste0("Fiber-seq_", id, "_sm.pdf"), 2, 3)
  }
}

main()
