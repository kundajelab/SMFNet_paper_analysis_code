#Run from the repository root: Rscript figure_code/Fig5_E8_F1_analysis/plot_disrupted_examples.R
dir.create("output/Fig5_E8_F1_analysis/plots", recursive=TRUE, showWarnings=FALSE)
PLOT_DIR <- "output/Fig5_E8_F1_analysis/plots/"

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(forcats)
filter <- dplyr::filter
library(hdf5r)
library(egg)
library(gridExtra)

theme_set(theme_classic() +
            theme(
              legend.text = element_text(size = 6),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
              axis.line = element_blank(),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 6),
              title = element_text(size = 6),
              strip.text = element_text(size = 6),
              strip.background = element_blank()
            ))

BASE <- "output/Fig5_E8_F1_analysis"

H5_FILE <- paste0(BASE, "/F1_negative_delta_expanded_400bp_windows.h5")
TFBS_FILE <- "input_files/Fig5_E8_F1_analysis/F1_annotated_motif_hits.tsv"
PLOT_DIR <- paste0(BASE, "/plots/variant_accessible_sorted")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

MOTIFS <- c("SP-KLF", "OCT4::SOX2", "OCT4::SOX2-alt")
CONDITIONS <- c("actual_b6", "actual_spret")
CONDITION_LABELS <- c(actual_b6 = "B6", actual_spret = "SPRET")

READ_COLORS <- c("0" = "#DEDEDE", "1" = "#000000")
TFBS_COLORS <- c("No variant" = "#B5B5B5", "Variant" = "#C65D65")
#TFBS annotations
tfbs_file <- read.delim(TFBS_FILE, stringsAsFactors = FALSE) %>%
  mutate(
    motif_start = as.integer(motif_start),
    motif_end = as.integer(motif_end),
    variant_status = if_else(
      !is.na(motif_variant) &
        !tolower(trimws(motif_variant)) %in% c("none", "", "na"),
      "Variant", "No variant"
    )
  )

tfbs_by_region <- split(tfbs_file, tfbs_file$region_id)

get_window_tfbs <- function(region, window_start, window_end) {
  hits <- tfbs_by_region[[region]]
  if (is.null(hits)) return(data.frame())
  
  hits <- hits %>%
    filter(motif_end > window_start, motif_start < window_end) %>%
    distinct(motif_start, motif_end, motif, motif_variant, source_allele, .keep_all = TRUE) %>%
    arrange(motif_start, motif_end)
  
  if (!nrow(hits)) return(hits)
  
  hits %>%
    mutate(
      xmin = pmax(-200.5, motif_start - window_start - 200 - 0.5),
      xmax = pmin(199.5, motif_end - window_start - 200 - 0.5),
      center = (xmin + xmax) / 2,
      lane = (row_number() - 1L) %% 4L,
      label_y = 0.095 + lane * 0.055
    )
}
#Helpers
read_matrix <- function(g, dataset) {
  mat <- g[[dataset]][, ]
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = 1)
  
  if (ncol(mat) != 400) {
    if (nrow(mat) == 400) mat <- t(mat)
    else stop(dataset, ": expected reads x 400 matrix")
  }
  mat
}

#Width of the continuous accessible region containing x = 0.
accessible_width <- function(row, center_idx) {
  center <- center_idx[which.min(abs(center_idx - mean(center_idx)))]
  
  left_nuc <- which(row[seq_len(center)] == 1)
  left <- if (length(left_nuc)) max(left_nuc) + 1L else 1L
  
  right_nuc <- which(row[center:length(row)] == 1)
  right <- if (length(right_nuc)) center + min(right_nuc) - 2L else length(row)
  
  max(0L, right - left + 1L)
}

#Approximate similarity ordering for inaccessible reads.
#Every 4th bp is used to make hierarchical clustering quicker.
order_inaccessible <- function(mat) {
  if (nrow(mat) <= 2) return(seq_len(nrow(mat)))
  
  sig <- mat[, seq(1, ncol(mat), by = 4), drop = FALSE]
  d <- dist(sig, method = "manhattan")
  hclust(d, method = "average")$order
}

#Accessible = entire central 15 bp is non-nucleosomal.
order_reads <- function(mat, x) {
  center_idx <- which(abs(x) <= 7.5)
  
  accessible <- rowSums(mat[, center_idx, drop = FALSE] == 1) == 0
  acc_idx <- which(accessible)
  inac_idx <- which(!accessible)
  
  #Accessible: largest accessible region first.
  if (length(acc_idx)) {
    widths <- vapply(
      acc_idx,
      function(i) accessible_width(mat[i, ], center_idx),
      numeric(1)
    )
    acc_idx <- acc_idx[order(widths, decreasing = TRUE)]
  }
  
  #Inaccessible: group approximately by profile similarity.
  if (length(inac_idx) > 1) {
    inac_order <- order_inaccessible(mat[inac_idx, , drop = FALSE])
    inac_idx <- inac_idx[inac_order]
  }
  
  list(
    order = c(acc_idx, inac_idx),
    accessible = accessible,
    fraction = mean(accessible),
    n_accessible = length(acc_idx),
    n_total = nrow(mat)
  )
}
#Average profile
make_average_df <- function(matrices, x, source_allele) {
  source_suffix <- if (source_allele == "B6") "b6" else "spret"
  
  bind_rows(lapply(CONDITIONS, function(cond) {
    allele <- if (grepl("_b6$", cond)) "b6" else "spret"
    
    data.frame(
      condition = cond,
      pos = x,
      signal = colMeans(matrices[[cond]], na.rm = TRUE),
      allele_role = ifelse(
        allele == source_suffix,
        "Source allele",
        "Opposite allele"
      )
    )
  }))
}

plot_average <- function(df, hits, region) {
  p <- ggplot(df, aes(pos, signal, color = allele_role)) +
    geom_line(linewidth = 0.65) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "#777777",
      linewidth = 0.3
    ) +
    scale_color_manual(values = c(
      "Source allele" = "#111111",
      "Opposite allele" = "#C64B55"
    )) +
    scale_x_continuous(
      limits = c(-200.5, 199.5),
      breaks = c(-200, -100, 0, 100, 199),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      x = "Position relative to motif center (bp)",
      y = "Mean nucleosome occupancy",
      color = NULL,
      title = region
    ) +
    theme(legend.position = "bottom")
  
  if (nrow(hits)) {
    p <- p +
      geom_rect(
        data = hits,
        aes(
          xmin = xmin, xmax = xmax,
          ymin = 0.025, ymax = 0.065,
          fill = variant_status
        ),
        inherit.aes = FALSE,
        color = NA
      ) +
      geom_text(
        data = hits,
        aes(
          x = center, y = label_y,
          label = motif,
          color = variant_status
        ),
        inherit.aes = FALSE,
        size = 4/2.8,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = TFBS_COLORS, name = "TFBS")
  }
  
  p
}
#Single-molecule plot
make_single_molecule_df <- function(matrices, x) {
  molecule_list <- list()
  stats_list <- list()
  
  for (cond in CONDITIONS) {
    mat <- matrices[[cond]]
    info <- order_reads(mat, x)
    
    mat <- mat[info$order, , drop = FALSE]
    n <- nrow(mat)
    
    #First row is at the TOP.
    y <- 1 - ((seq_len(n) - 0.5) / n)
    
    molecule_list[[cond]] <- data.frame(
      condition = cond,
      pos = rep(x, times = n),
      y = rep(y, each = length(x)),
      state = as.vector(t(mat))
    )
    
    stats_list[[cond]] <- data.frame(
      condition = cond,
      fraction = info$fraction,
      n_accessible = info$n_accessible,
      n_total = info$n_total,
      boundary_y = 1 - info$fraction
    )
  }
  
  list(
    molecules = bind_rows(molecule_list) %>%
      mutate(
        condition = factor(condition, CONDITIONS, CONDITION_LABELS),
        state = factor(state, levels = c(0, 1))
      ),
    stats = bind_rows(stats_list) %>%
      mutate(
        condition = factor(condition, CONDITIONS, CONDITION_LABELS),
        label = paste0(
          round(100 * fraction, 1),
          "% accessible"
        )
      )
  )
}

plot_single_molecule <- function(sm) {
  ggplot(sm$molecules, aes(pos, y, fill = state)) +
    geom_raster(interpolate = FALSE) +
    
    #Boundary between accessible and inaccessible reads.
    geom_hline(
      data = sm$stats,
      aes(yintercept = boundary_y),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.3,
      color = "#B65B60"
    ) +
    
    #Accessible fraction.
    geom_text(
      data = sm$stats,
      aes(x = -192, y = 0.97, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.5
    ) +
    
    facet_wrap(~condition, nrow = 1) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#B65B60",
      linewidth = 0.25
    ) +
    
    scale_fill_manual(values = READ_COLORS, drop = FALSE) +
    
    scale_x_continuous(
      limits = c(-200.5, 199.5),
      breaks = c(-200, 0, 199),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      expand = c(0, 0)
    ) +
    
    labs(
      x = "Relative position (bp)",
      y = "Single molecules"
    ) +
    
    theme(legend.position = "none")
}
#Process windows
h5 <- H5File$new(H5_FILE, mode = "r")
cw <- h5[["windows"]]

#Filter directly from motif attributes in the H5.
window_ids <- names(cw)[
  vapply(
    names(cw),
    function(id) {
      g <- cw[[id]]
      motif <- as.character(h5attr(g, "motif"))
      motif %in% MOTIFS
    },
    logical(1)
  )
]

#The two experimental examples in Fig5i-j and E8h-i.
window_ids <- window_ids[vapply(window_ids, function(id) {
  g <- cw[[id]]
  region <- as.character(h5attr(g, 'region_id'))
  start <- as.integer(h5attr(g, 'window_start'))
  (region == 'chr1:46904135-46906183' && start == 46905014L) ||
    (region == 'chr7:27621268-27623316' && start == 27622166L)
}, logical(1))]

cat("Selected windows:", length(window_ids), "\n")

for (i in seq_along(window_ids)) {
  window_id <- window_ids[i]
  g <- cw[[window_id]]
  
  region <- as.character(h5attr(g, "region_id"))
  motif <- as.character(h5attr(g, "motif"))
  window_start <- as.integer(h5attr(g, "window_start"))
  window_end <- as.integer(h5attr(g, "window_end"))
  
  source_allele <- toupper(as.character(h5attr(g, "source_allele")))
  if (source_allele == "BL6") source_allele <- "B6"
  
  cat(sprintf(
    "\rProcessing %d/%d: %s | %s",
    i, length(window_ids), region, motif
  ))
  flush.console()
  
  genomic_positions <- as.integer(g[["genomic_positions"]][])
  x <- genomic_positions - window_start - 200L
  stopifnot(length(x) == 400)
  
  matrices <- setNames(
    lapply(CONDITIONS, function(cond) read_matrix(g, cond)),
    CONDITIONS
  )
  
  hits <- get_window_tfbs(region, window_start, window_end)
  
  avg_df <- make_average_df(matrices, x, source_allele)
  avg_plot <- plot_average(avg_df, hits, region)
  
  sm <- make_single_molecule_df(matrices, x)
  sm_plot <- plot_single_molecule(sm)
  
  safe <- gsub(
    "[:\\\\/|-]",
    "_",
    paste(region, motif, window_start, window_end, sep = "_")
  )
  
  pdf(
    file.path(PLOT_DIR, paste0(safe, "_avg.pdf")),
    width = 20/2.54,
    height = 20/2.54
  )
  grid.arrange(set_panel_size(avg_plot, height = unit(2, "cm"), width = unit(3, "cm")))
  dev.off()
  
  pdf(
    file.path(PLOT_DIR, paste0(safe, "_single_molecules.pdf")),
    width = 20/2.54,
    height = 20/2.54
  )
  grid.arrange(set_panel_size(sm_plot, height = unit(2, "cm"), width = unit(2, "cm")))
  dev.off()
}

cat("\n")
h5$close_all()
