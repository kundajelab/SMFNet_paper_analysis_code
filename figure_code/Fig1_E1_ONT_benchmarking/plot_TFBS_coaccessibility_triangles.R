#Run from the repository root: Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_TFBS_coaccessibility_triangles.R
#Run prepare_TFBS_pairs.py and compute_TFBS_coaccessibility.py first.
#Plots triangles above experiment/model single-molecule states for each eligible region.
library(ggplot2)
library(dplyr)
library(tidyr)
library(egg)
library(gridExtra)
library(scales)
library(hdf5r)

input_dir <- file.path("input_files", "Fig1_E3_Codependency_benchmarking")
analysis_dir <- file.path("output", "Fig1_E3_Codependency_benchmarking", "TFBS")
output_dir <- file.path(analysis_dir, "triangles")
experiment_states <- file.path(input_dir, "mESC_SMF-ONT_testRegions_Experiment_states.h5")
model_states <- file.path(input_dir, "mESC_SMF-ONT_testRegions_SMFNet_states.h5")
co_ac <- read.delim(file.path(analysis_dir, "TFBS_coaccessibility.tsv"))
correction <- read.delim(file.path(input_dir, "coaccessibility_background_bias.txt"))
bias <- correction$bias[match(c("mi_act", "mi_pred"), correction$type)]
stopifnot(length(bias) == 2, all(is.finite(bias)))
co_ac <- co_ac %>% mutate(mi_act = mi_act - bias[1], mi_pred = mi_pred - bias[2])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), plot.title = element_text(size = 6), axis.text = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6), strip.text = element_text(size = 6),
  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
))

#Original midpoint interpolation and split clustering from plot_smf_tfbs_coac.
extend_to_full <- function(mat, smf_pos, full_len = 512) {
  
  mat <- as.matrix(mat)
  smf_pos <- as.numeric(smf_pos)
  
  n <- nrow(mat)
  out <- matrix(0, nrow = n, ncol = full_len)
  
  for (r in seq_len(n)) {
    
    row <- mat[r, ]
    
    #Treat only NA as missing (keeps -1)
    valid <- !is.na(row)
    idx <- which(valid)
    
    if (length(idx) == 0)
      next
    
    pos <- smf_pos[idx]
    val <- row[idx]
    
    keep <- pos >= 1 & pos <= full_len
    
    mapped_pos <- as.integer(pos[keep])
    val <- val[keep]
    
    if (length(mapped_pos) == 0)
      next
    
    ord <- order(mapped_pos)
    mapped_pos <- mapped_pos[ord]
    val <- val[ord]
    
    #left fill
    out[r, seq_len(mapped_pos[1])] <- val[1]
    
    #interpolate
    for (i in seq_len(length(mapped_pos) - 1)) {
      
      p1 <- mapped_pos[i]
      p2 <- mapped_pos[i + 1]
      
      if (p2 <= p1)
        next
      
      if (val[i] == val[i + 1]) {
        
        out[r, p1:(p2 - 1)] <- val[i]
        
      } else {
        
        mid <- floor((p1 + p2) / 2)
        
        out[r, p1:(mid - 1)] <- val[i]
        out[r, mid:(p2 - 1)] <- val[i + 1]
      }
    }
    
    #right fill
    out[r, mapped_pos[length(mapped_pos)]:full_len] <-
      val[length(val)]
  }
  
  out
}

cluster_split <- function(states_pred, states_act, split_col = 200, k = 15){
    
    cluster_group <- function(pred, act){
      
      both <- rbind(pred, act)
      
      if(nrow(both) <= 1){
        return(list(
          pred = seq_len(nrow(pred)),
          act = seq_len(nrow(act))
        ))
      }
      
      hc <- hclust(
        dist(both),
        method = "ward.D2"
      )
      
      k_use <- min(k, nrow(both))
      
      cl <- cutree(hc, k = k_use)
      
      cl_mean <- tapply(
        rowMeans(both),
        cl,
        mean
      )
      
      cl_order <- as.integer(names(sort(cl_mean, decreasing = FALSE)))
      
      dend_order <- hc$order
      
      new_order <- integer(0)
      
      for(cc in cl_order){
        
        rows <- which(cl == cc)
        
        #preserve dendrogram ordering within cluster
        rows <- dend_order[dend_order %in% rows]
        
        new_order <- c(new_order, rows)
      }
      
      list(
        pred = new_order[new_order <= nrow(pred)],
        act = new_order[new_order > nrow(pred)] - nrow(pred)
      )
    }
    
    pred0 <- which(states_pred[, split_col] == 0)
    act0 <- which(states_act[, split_col] == 0)
    
    pred1 <- which(states_pred[, split_col] == 1)
    act1 <- which(states_act[, split_col] == 1)
    
    ord0 <- cluster_group(
      states_pred[pred0, , drop = FALSE],
      states_act[act0, , drop = FALSE]
    )
    
    ord1 <- cluster_group(
      states_pred[pred1, , drop = FALSE],
      states_act[act1, , drop = FALSE]
    )
    
    list(
      pred = c(pred0[ord0$pred], pred1[ord1$pred]),
      act = c(act0[ord0$act], act1[ord1$act])
    )
  }
  


#Use the original state interpolation and read filter over the full triangle window.
read_region_states <- function(path, region) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  positions <- as.numeric(h5[["smf_pos"]][[region]][])
  states <- t(h5[["smf_states"]][[region]]$read(drop = FALSE))
  stopifnot(ncol(states) == length(positions))
  keep <- positions >= 0 & positions <= 1024
  positions <- positions[keep]
  states <- states[, keep, drop = FALSE]
  states <- states[rowMeans(states) != -1, , drop = FALSE]
  states <- extend_to_full(states, positions, full_len = 1024)
  #States 0/1/2 are accessible; state 3 is non-accessible; -1 remains missing.
  states[states %in% c(0, 1, 2)] <- 0
  states[states == 3] <- 1
  states[between(rowMeans(states), 0.05, 0.95), , drop = FALSE]
}

molecule_plot <- function(region, motifs) {
  act <- read_region_states(experiment_states, region)
  pred <- read_region_states(model_states, region)
  #Preserve the original midpoint split, five clusters, and within-cluster ordering.
  ordering <- cluster_split(pred, act, split_col = 512, k = 5)
  act <- act[ordering$act, , drop = FALSE]
  pred <- pred[ordering$pred, , drop = FALSE]
  make_table <- function(matrix, type) {
    data.frame(pos = rep(1:1024, times = nrow(matrix)),
               molecule = rep(seq_len(nrow(matrix)), each = 1024),
               state = factor(as.vector(t(matrix)), levels = c(-1, 0, 1)), type = rep(type, nrow(matrix) * 1024))
  }
  data <- bind_rows(make_table(act, "mi_act"), make_table(pred, "mi_pred"))
  ggplot(data, aes(pos, -molecule, fill = state)) +
    geom_raster() + facet_wrap(~type, scales = "free_y", labeller = as_labeller(c(mi_act = "Experiment", mi_pred = "SMFNet"))) +
    geom_vline(data = motifs, aes(xintercept = pos), color = "#888888",
               linewidth = 0.15, linetype = "dashed") +
    scale_x_continuous(breaks = c(0, 512, 1024), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1024)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("-1" = "#FFFFFF", "0" = "#DEDEDE", "1" = "#000000"),
                      labels = c("Missing", "Accessible", "Non-accessible"), drop = FALSE) +
    labs(x = "Relative position (bp)", y = "Single molecules", fill = "State") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

for (region in unique(co_ac$region)) {
  region_data <- co_ac[co_ac$region == region, ]
  motifs <- bind_rows(
    transmute(region_data, motif = motif1, pos = pos1),
    transmute(region_data, motif = motif2, pos = pos2)
  ) %>% distinct()
  #Retain the original requirement of more than five annotated TFBS per plotted region.
  if (nrow(motifs) <= 5) next
  data <- region_data %>%
    pivot_longer(c(mi_act, mi_pred), names_to = "type", values_to = "mi") %>%
    mutate(mid_pos = (pos1 + pos2) / 2)
  plot <- ggplot(data, aes(mid_pos, distance, fill = mi)) +
    facet_wrap(~type, labeller = as_labeller(c(mi_act = "Experiment", mi_pred = "SMFNet"))) +
    geom_rect(data = motifs, inherit.aes = FALSE,
      aes(xmin = pos - 10, xmax = pos + 10, ymin = 10, ymax = 50), fill = "#DEDEDE") +
    geom_text(data = motifs, inherit.aes = FALSE, aes(pos, 30, label = motif),
      size = 4 / 2.8, vjust = 1.2) +
    scale_x_continuous(limits = c(0, 1024), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1024)) +
    geom_point(shape = 21, size = 1.5) +
    scale_fill_gradient2(low = "#008837", mid = "white", high = "#7b3294",
      trans = pseudo_log_trans(sigma = 1e-3), limits = c(-1, 1),
      breaks = c(1, 0.1, 0.01, 0, -0.01, -0.1, -1)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
      axis.line.x = element_line(color = "black"), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.title.y = element_blank())
  path <- file.path(output_dir, paste0(gsub("[:-]", "_", region), "_TFBS_coaccessibility.pdf"))
  states_plot <- molecule_plot(region, motifs)
  triangle_panel <- set_panel_size(plot, height = unit(2, "cm"), width = unit(3, "cm"))
  molecule_panel <- set_panel_size(states_plot, height = unit(3, "cm"), width = unit(3, "cm"))
  #Align the experiment/model columns and the genomic x axes across both rows.
  widths <- grid::unit.pmax(triangle_panel$widths, molecule_panel$widths)
  triangle_panel$widths <- widths
  molecule_panel$widths <- widths
  pdf(path, width = 20 / 2.54, height = 20 / 2.54)
  grid.arrange(triangle_panel, molecule_panel, ncol = 1)
  dev.off()
  message("Wrote: ", path)
}
