#Run from the repository root: Rscript figure_code/Fig5_E8_F1_analysis/plot_microstate_examples.R
dir.create("output/Fig5_E8_F1_analysis/plots", recursive=TRUE, showWarnings=FALSE)
PLOT_DIR <- "output/Fig5_E8_F1_analysis/plots/"

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(forcats)
filter <- dplyr::filter
library(egg)
library(gridExtra)
library(data.table)
library(viridis)
library(hdf5r)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6), 
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

PLOT_DIR <- "output/Fig5_E8_F1_analysis/plots/"
BASE <- "output/Fig5_E8_F1_analysis"
EXPANDED_H5 <- paste0(BASE, "/F1_expanded_400bp_windows.h5")
CLUSTER_H5 <- paste0(BASE, "/F1_all4_balanced_clusters.h5")
WINDOW_DF <- paste0(BASE, "/F1_ac_pred_act.tsv")
TFBS_FILE <- "input_files/Fig5_E8_F1_analysis/F1_annotated_motif_hits.tsv"
PLOT_DIR <- paste0(BASE, "/plots/f1_cluster_plots")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

CONDITIONS <- c("actual_b6", "actual_spret", "predicted_b6", "predicted_spret")
CONDITION_LABELS <- c(actual_b6 = "Actual B6", actual_spret = "Actual SPRET",
                      predicted_b6 = "Predicted B6", predicted_spret = "Predicted SPRET")
READ_COLORS <- c("0" = "#DEDEDE", "1" = "#000000", "-1" = "#C8C8C8")
TFBS_COLORS <- c("No variant" = "#B5B5B5", "Variant" = "#C65D65")
CLUSTER_LETTERS <- LETTERS[1:8]
#Load annotations
#Motif-centered windows, including the focal source allele.
window_df <- read.delim(WINDOW_DF, stringsAsFactors = FALSE) %>% 
  filter(actual_delta_accessible >= 0.1) %>% 
  filter(abs(actual_delta_accessible - predicted_delta_accessible) <= 0.05) %>% 
  filter(motif %in% c("OCT4::SOX2", "SP-KLF", "MAX::MYC", "NRF1", "BANP"))
stopifnot(all(c("h5_window_id", "source_allele") %in% names(window_df)))

#All TFBS, including motifs without variants.
window_df <- window_df %>% filter(region_id == 'chr6:124399174-124401222')

tfbs_file <- read.delim(TFBS_FILE, stringsAsFactors = FALSE) %>%
  mutate(motif_start = as.integer(motif_start), motif_end = as.integer(motif_end),
         variant_status = if_else(
           !is.na(motif_variant) &
             !tolower(trimws(motif_variant)) %in% c("none", "", "na"),
           "Variant", "No variant"
         ))
#Helpers
read_matrix <- function(group, dataset, n_labels) {
  mat <- group[[dataset]][, ]
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = 1)
  
  if (nrow(mat) != n_labels) {
    if (ncol(mat) == n_labels) mat <- t(mat)
    else stop(dataset, ": matrix dimensions do not match cluster labels")
  }
  
  stopifnot(ncol(mat) == 400)
  mat
}

get_source_allele <- function(window_id) {
  matches <- window_df %>% filter(h5_window_id == window_id)
  alleles <- unique(toupper(trimws(na.omit(matches$source_allele))))
  alleles[alleles == "BL6"] <- "B6"
  alleles <- unique(alleles)
  
  if (length(alleles) != 1 || !alleles %in% c("B6", "SPRET")) {
    stop(window_id, ": cannot uniquely determine source allele from window_df")
  }
  
  alleles
}

get_window_tfbs <- function(region, window_start, window_end) {
  hits <- tfbs_file %>%
    filter(region_id == region, motif_end > window_start, motif_start < window_end) %>%
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
#Average occupancy: all clusters combined
make_average_df <- function(matrices, x, source_allele) {
  source_suffix <- if (source_allele == "B6") "b6" else "spret"
  
  bind_rows(lapply(CONDITIONS, function(condition) {
    allele <- if (grepl("_b6$", condition)) "b6" else "spret"
    
    data.frame(
      condition = condition,
      pos = x,
      signal = colMeans(matrices[[condition]], na.rm = TRUE),
      allele_role = if (allele == source_suffix) "Source allele" else "Opposite allele",
      data_type = if (startsWith(condition, "actual")) "Actual" else "Predicted"
    )
  })) %>%
    mutate(
      allele_role = factor(allele_role, levels = c("Source allele", "Opposite allele")),
      data_type = factor(data_type, levels = c("Actual", "Predicted"))
    )
}

plot_average <- function(avg_df, hits, source_allele, region, window_start, window_end) {
  p <- ggplot(avg_df, aes(x = pos, y = signal, color = allele_role)) +
    facet_wrap(~data_type) +
    geom_line(linewidth = 0.65) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "#777777", linewidth = 0.3) +
    scale_color_manual(values = c("Source allele" = "#111111",
                                  "Opposite allele" = "#C64B55")) +
    scale_x_continuous(limits = c(-200.5, 199.5),
                       breaks = c(-200, -100, 0, 100, 199), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      x = "Position relative to motif center (bp)",
      y = "Mean occupancy",
      color = NULL,
      linetype = NULL,
      title = paste0(region, " | ", window_start, "_", window_end, " | REF: ", source_allele)
    ) +
    theme(legend.position = "bottom")
  
  #Motif boxes lie near y = 0.05, with names staggered above them.
  #Draw these after the curves so every annotated motif remains visible.
  if (nrow(hits)) {
    p <- p +
      geom_rect(
        data = hits,
        aes(xmin = xmin, xmax = xmax, ymin = 0.025, ymax = 0.065,
            fill = variant_status),
        inherit.aes = FALSE, color = NA, alpha = 0.95
      ) +
      geom_text(
        data = hits,
        aes(x = center, y = label_y, label = motif, color = variant_status),
        inherit.aes = FALSE, size = 4/2.8, check_overlap = TRUE, show.legend = FALSE
      ) +
      scale_fill_manual(values = TFBS_COLORS, name = "TFBS") +
      guides(color = guide_legend(order = 1),
             linetype = guide_legend(order = 2),
             fill = guide_legend(order = 3))
  }
  
  p
}
#Single-molecule data: all reads, no sampling
make_single_molecule_df <- function(matrices, labels, x, cluster_levels) {
  result <- list()
  fractions <- list()
  
  for (condition in CONDITIONS) {
    mat <- matrices[[condition]]
    cl <- labels[[condition]]
    total_reads <- length(cl)
    
    for (k in cluster_levels) {
      members <- which(cl == k)
      fraction <- length(members) / total_reads
      
      fractions[[length(fractions) + 1L]] <- data.frame(
        condition = condition, cluster = CLUSTER_LETTERS[k + 1L],
        n_cluster = length(members), total_reads = total_reads, fraction = fraction
      )
      
      #Every facet gets exactly total_reads rows. Only the bottom
      #n_cluster rows are occupied; the remainder is transparent/white.
      raster_mat <- matrix(NA_integer_, nrow = total_reads, ncol = length(x))
      
      if (length(members)) {
        subset_mat <- mat[members, , drop = FALSE]
        subset_mat <- subset_mat[order(rowMeans(subset_mat, na.rm = TRUE)), , drop = FALSE]
        raster_mat[seq_len(nrow(subset_mat)), ] <- subset_mat
      }
      
      result[[length(result) + 1L]] <- data.frame(
        condition = condition, cluster = CLUSTER_LETTERS[k + 1L],
        pos = rep(x, times = total_reads),
        y = rep((seq_len(total_reads) - 0.5) / total_reads, each = length(x)),
        state = as.vector(t(raster_mat))
      )
    }
  }
  
  molecule_df <- bind_rows(result) %>%
    mutate(condition = factor(condition, levels = CONDITIONS, labels = CONDITION_LABELS),
           cluster = factor(cluster, levels = CLUSTER_LETTERS[cluster_levels + 1L]),
           state = factor(state, levels = c(0, 1, -1)))
  
  fraction_df <- bind_rows(fractions) %>%
    mutate(condition = factor(condition, levels = CONDITIONS, labels = CONDITION_LABELS),
           cluster = factor(cluster, levels = CLUSTER_LETTERS[cluster_levels + 1L]))
  
  list(molecules = molecule_df, fractions = fraction_df)
}

plot_single_molecule <- function(molecule_df, fraction_df, source_allele) {
  fraction_labels <- fraction_df %>%
    mutate(label = paste0(round(100 * fraction, 1), "% (", n_cluster, "/", total_reads, ")"))
  
  ggplot(molecule_df, aes(x = pos, y = y, fill = state)) +
    geom_raster(interpolate = FALSE) +
    geom_text(data = fraction_labels, aes(x = -192, y = 0.45, label = label),
              inherit.aes = FALSE, hjust = 0, size = 2.4) +
    facet_grid(cluster ~ condition, drop = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "#B65B60", linewidth = 0.25) +
    scale_fill_manual(values = READ_COLORS,
                      na.value = "white", drop = FALSE) +
    scale_x_continuous(limits = c(-200.5, 199.5),
                       breaks = c(-200, 0, 199), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.5),
                       expand = c(0, 0)) +
    labs(x = "Relative position (bp)", y = "Single molecules", 
         title = paste0("Single molecules by cluster | Source allele: ", source_allele)) +
    theme(legend.position = "none")
}
#Open HDF5 files and organize windows by region
expanded_h5 <- H5File$new(EXPANDED_H5, mode = "r")
cluster_h5 <- H5File$new(CLUSTER_H5, mode = "r")

expanded_windows <- expanded_h5[["windows"]]
cluster_windows <- cluster_h5[["windows"]]

window_ids <- intersect(names(expanded_windows), names(cluster_windows))
window_ids <- intersect(window_ids, window_df$h5_window_id)

window_regions <- vapply(
  window_ids,
  function(id) as.character(hdf5r::h5attr(expanded_windows[[id]], "region_id")),
  character(1)
)

windows_by_region <- split(window_ids, window_regions)

cat("Clustered windows:", length(window_ids),
    "| Regions:", length(windows_by_region), "\n")
#Plot every clustered window
for (i in seq_along(windows_by_region)) {
  region <- names(windows_by_region)[i]
  region_dir <- file.path(PLOT_DIR, gsub("[:-]", "_", region))
  dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\rProcessing region %d/%d: %s", i, length(windows_by_region), region))
  flush.console()
  
  for (window_id in windows_by_region[[region]]) {
    g_exp <- expanded_windows[[window_id]]
    g_cl <- cluster_windows[[window_id]]
    
    motif <- window_df[window_df$region == region,]$motif
    safe_region <- gsub("[:\\-]", "_", region)
    window_start <- as.integer(hdf5r::h5attr(g_exp, "window_start"))
    window_end <- as.integer(hdf5r::h5attr(g_exp, "window_end"))
    motif_start <- as.integer(hdf5r::h5attr(g_exp, "motif_start"))
    motif_end <- as.integer(hdf5r::h5attr(g_exp, "motif_end"))
    
    genomic_positions <- as.integer(g_exp[["genomic_positions"]][])
    x <- genomic_positions - window_start - 200L
    stopifnot(length(x) == 400)
    
    source_allele <- get_source_allele(window_id)
    
    labels <- setNames(
      lapply(CONDITIONS, function(cond)
        as.integer(g_cl[[paste0("labels_", cond)]][])),
      CONDITIONS
    )
    
    matrices <- setNames(
      lapply(CONDITIONS, function(cond)
        read_matrix(g_exp, cond, length(labels[[cond]]))),
      CONDITIONS
    )
    
    cluster_levels <- sort(unique(unlist(labels, use.names = FALSE)))
    cluster_levels <- cluster_levels[!is.na(cluster_levels) & cluster_levels >= 0]
    
    if (!length(cluster_levels)) next
    if (any(cluster_levels >= 8)) stop(window_id, ": cluster labels exceed A-H")
    
    hits <- get_window_tfbs(region, window_start, window_end)
    avg_df <- make_average_df(matrices, x, source_allele)
    sm_data <- make_single_molecule_df(matrices, labels, x, cluster_levels)
    
    avg_plot <- plot_average(avg_df, hits, source_allele, region, window_start, window_end)
    sm_plot <- plot_single_molecule(
      sm_data$molecules,
      sm_data$fractions,
      source_allele
    )
    
    
    pdf(paste0(region_dir, "/", safe_region, "_", motif, "_avg.pdf"),width=20/2.54,height=20/2.54)
    grid.arrange(set_panel_size(avg_plot, height=unit(2,"cm"),width=unit(3,"cm"))); dev.off()
    
    
    pdf(paste0(region_dir, "/", safe_region, "_", motif, "_sm_cluster.pdf"),width=20/2.54,height=20/2.54)
    grid.arrange(set_panel_size(sm_plot, height=unit(0.75,"cm"),width=unit(2.5,"cm"))); dev.off()

  }
}
#Close files and save plot log
expanded_h5$close_all()
cluster_h5$close_all()

#Plot cluster changes
cluster_delta_dir <- paste0(BASE, "/plots/example_cluster_delta/")
dir.create(cluster_delta_dir, recursive=TRUE, showWarnings=FALSE)

cluster_h5 <- H5File$new(CLUSTER_H5, mode = "r")
cw <- cluster_h5[["windows"]]

cluster_frac <- bind_rows(lapply(names(cw), function(window_id) {
  g <- cw[[window_id]]
  region <- as.character(h5attr(g, "region_id"))
  
  bind_rows(lapply(CONDITIONS, function(cond) {
    lab <- as.integer(g[[paste0("labels_", cond)]][])
    as.data.frame(table(cluster = lab)) %>%
      mutate(
        window_id = window_id,
        region_id = region,
        condition = cond,
        fraction = Freq / sum(Freq)
      )
  }))
}))

cluster_h5$close_all()

cluster_frac_fil <- cluster_frac %>%
  select(region_id, window_id, condition, cluster, count = Freq, fraction) %>% 
  filter(region_id %in% window_df$region_id) %>% 
  mutate(cluster = LETTERS[as.integer(as.character(cluster)) + 1])

for (r in unique(cluster_frac_fil$region_id)) {
  safe_region <- gsub("[:\\-]", "_", r)
  
  loop_df <- cluster_frac_fil %>% 
    filter(region_id == r) %>% 
    select(-count) %>% 
    pivot_wider(names_from = condition, values_from = fraction) %>% 
    mutate(act_delta = actual_spret - actual_b6, pred_delta = predicted_spret - predicted_b6)
  
  
  delta_df <- loop_df %>% 
    select(region_id, cluster, act_delta, pred_delta) %>% 
    pivot_longer(-c(region_id, cluster), names_to = "type", values_to = "delta")
  
  delta_plot <- delta_df %>% 
    ggplot() +
    facet_wrap(~type, ncol = 1) +
    geom_bar(aes(y = cluster, x = delta, fill = cluster), stat = "identity") +
    scale_x_continuous(limits = c(-0.2, 0.2)) +
    scale_y_discrete(limits = rev) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    scale_fill_manual(values = grey.colors(8, start = 0.2, end = 0.8)) +
    theme(axis.ticks.y = element_blank(), legend.position = "none")
  
  pdf(paste0(cluster_delta_dir, safe_region, "_cluster_perc.pdf"),width=20/2.54,height=20/2.54)
  grid.arrange(set_panel_size(delta_plot, height=unit(2,"cm"),width=unit(1.5,"cm"))); dev.off()
  

  }
