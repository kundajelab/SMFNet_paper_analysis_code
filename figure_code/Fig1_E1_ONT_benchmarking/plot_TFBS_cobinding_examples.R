#Run from the repository root: Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_TFBS_cobinding_examples.R
#Run prepare_TFBS_pairs.py and compute_TFBS_cobinding.py first.
#The original mean profiles use states: 0/1 -> 0, 2/3 -> 1, excluding -1.
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(egg)
library(gridExtra)
library(hdf5r)

input_dir <- file.path("input_files", "Fig1_E3_Codependency_benchmarking")
analysis_dir <- file.path("output", "Fig1_E3_Codependency_benchmarking", "TFBS", "cobinding")
output_dir <- file.path(analysis_dir, "examples")
pred_h5 <- file.path(input_dir, "mESC_SMF-ONT_testRegions_SMFNet_states.h5")
act_h5 <- file.path(input_dir, "mESC_SMF-ONT_testRegions_Experiment_states.h5")
co_binding <- read.delim(file.path(analysis_dir, "TFBS_cobinding.tsv"))
#Retain the original example selection and the original complete-case plotting cohort.
fil_tfbs <- co_binding %>% filter(as.logical(plot_eligible)) %>%
  mutate(diff_odds = abs(log_odds_act - log_odds_pred)) %>%
  filter(log_odds_act >= 2, log_odds_pred >= 2, diff_odds <= 2, distance >= 30)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(fil_tfbs, file.path(output_dir, "selected_examples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(), axis.text = element_text(size = 6), plot.title = element_text(size = 6),
  title = element_text(size = 6), axis.title = element_text(size = 6), strip.text = element_text(size = 6),
  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")
))

extend_to_full <- function(mat, smf_pos, full_len = 1024) {

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
    out[r, mapped_pos[length(mapped_pos)]:full_len] <- val[length(val)]
  }

  out
}

cluster_split <- function(states_pred, states_act, pos1, pos2) {

  #desired ordering
  pair_levels <- c("2_2", "2_0", "2_3", "0_2", "3_2", "0_0", "0_3", "3_0", "3_3")

  #pair frequencies
  pair_pred <- paste0(states_pred[, as.character(pos1)], "_", states_pred[, as.character(pos2)])

  pair_act <- paste0(states_act[, as.character(pos1)], "_", states_act[, as.character(pos2)])

  state_fraction_pred <- prop.table(table(factor(pair_pred, levels = pair_levels)))

  state_fraction_act <- prop.table(table(factor(pair_act, levels = pair_levels)))

  #combine for clustering
  both <- rbind(states_pred, states_act)

  pair <- c(pair_pred, pair_act)

  if (nrow(both) <= 1) {
    return(list(pred = seq_len(nrow(states_pred)),
      act = seq_len(nrow(states_act)), state_fraction_pred = state_fraction_pred,
      state_fraction_act = state_fraction_act
    ))
  }

  pair <- factor(pair, levels = pair_levels, ordered = TRUE)

  new_order <- integer(0)

  for (lev in levels(pair)) {

    rows <- which(pair == lev)

    if (length(rows) == 0)
      next

    if (length(rows) == 1) {
      new_order <- c(new_order, rows)
      next
    }

    hc <- hclust(dist(both[rows, , drop = FALSE]), method = "ward.D2")

    new_order <- c(new_order, rows[hc$order])
  }

  list(pred = new_order[new_order <= nrow(states_pred)],
    act = new_order[new_order > nrow(states_pred)] - nrow(states_pred),
    state_fraction_pred = state_fraction_pred, state_fraction_act = state_fraction_act)
}

main <- function() {
  h5_pred_states <- H5File$new(pred_h5, mode = "r")
  on.exit(h5_pred_states$close_all(), add = TRUE)
  h5_act_states <- H5File$new(act_h5, mode = "r")
  on.exit(h5_act_states$close_all(), add = TRUE)
  for (i in seq_len(nrow(fil_tfbs))){
    region <- fil_tfbs[i,]$region
    pos1 <- fil_tfbs[i,]$pos1
    pos2 <- fil_tfbs[i,]$pos2
    motif1 <- fil_tfbs[i,]$motif1
    motif2 <- fil_tfbs[i,]$motif2
    mi_pred <- round(fil_tfbs[i,]$log_odds_pred, 2)
    mi_act <-  round(fil_tfbs[i,]$log_odds_act, 2)

    print(paste0("Processing [", i, "/", nrow(fil_tfbs), "] ", region))

    safe_region <- gsub("[:-]", "_", region)
    stem <- paste(safe_region, gsub("[^A-Za-z0-9_.-]", "_", motif1), pos1,
                  gsub("[^A-Za-z0-9_.-]", "_", motif2), pos2, sep = "_")
    window_start <- min(pos1, pos2) - 200
    window_end <- max(pos1, pos2) + 200

    #PLOT sm data

    states_pred <- t(h5_pred_states[["smf_states"]][[region]][,])
    states_act <- t(h5_act_states[["smf_states"]][[region]][,])

    pos_pred <- h5_pred_states[["smf_pos"]][[region]][]
    pos_act <- h5_act_states[["smf_pos"]][[region]][]

    colnames(states_pred) <- pos_pred
    colnames(states_act) <- pos_act

    keep_pred <- between(as.numeric(colnames(states_pred)), 0, 1024)
    keep_act  <- between(as.numeric(colnames(states_act)), 0, 1024)

    states_pred <- states_pred[, keep_pred, drop = FALSE]
    states_act  <- states_act[, keep_act, drop = FALSE]
    pos_pred <- pos_pred[keep_pred]
    pos_act <- pos_act[keep_act]

    #CALCULATE MEAN AVG FIRST
    binary_pred <- states_pred
    binary_pred[binary_pred %in% c(0,1)] <- 0
    binary_pred[binary_pred %in% c(2,3)] <- 1

    binary_act <- states_act
    binary_act[binary_act %in% c(0,1)] <- 0
    binary_act[binary_act %in% c(2,3)] <- 1

    mean_act <- apply(binary_act, 2, function(x){

      x <- x[x != -1]

      if(length(x) == 0)
        return(NA_real_)

      mean(x)

    })
    mean_pred <- apply(binary_pred, 2, function(x){

      x <- x[x != -1]

      if(length(x) == 0)
        return(NA_real_)

      mean(x)

    })

    act_avg_df <- data.frame(pos = names(mean_act), signal = mean_act, type = "act")

    pred_avg_df <- data.frame(pos = names(mean_pred), signal = mean_pred, type = "pred")

    avg_df <- rbind(pred_avg_df, act_avg_df)

    p_avg <- avg_df %>%
      filter(between(as.numeric(pos), window_start, window_end)) %>%
      mutate(pos = as.numeric(pos) - window_start +1) %>%
      ggplot(aes(x = pos, y = signal, color = type)) +
      geom_line() +
      geom_point(size = 0.5) +
      scale_y_continuous(limits = c(0,1), name = "Mean occupancy (states)") +
      scale_x_continuous(name = "Relative position (bp)", limits = c(1, window_end - window_start +1)) +
      scale_color_manual(values = c("#BCBCBC", "#636cb2")) +
      annotate("rect", xmin = pos1 - window_start - 7.5, xmax = pos1 - window_start + 7.5, ymin = 0, ymax = 0.1, fill = "#DEDEDE") +
      annotate("rect", xmin = pos2 - window_start - 7.5, xmax = pos2 - window_start + 7.5, ymin = 0, ymax = 0.1, fill = "#DEDEDE") +
      annotate("text", x = pos1 - window_start, y = 0.05, size = 5/2.8, label = motif1) +
      annotate("text", x = pos2 - window_start, y = 0.05, size = 5/2.8, label = motif2)

    pdf(file.path(output_dir, paste0(stem, "_average.pdf")), width = 20/2.54, height = 20/2.54)
    fix <- set_panel_size(p_avg, height = unit(2, "cm"), width = unit(3, "cm"))
    grid.arrange(fix)
    dev.off()

    #Now single molecule

    states_pred_full <- extend_to_full(states_pred, pos_pred)
    states_act_full <- extend_to_full(states_act, pos_act)
    colnames(states_pred_full) <- 1:1024
    colnames(states_act_full) <- 1:1024

    fil_pos_pred <- between(as.numeric(colnames(states_pred_full)), window_start, window_end)
    fil_pos_act <- between(as.numeric(colnames(states_act_full)), window_start, window_end)

    states_pred_fil <- states_pred_full[,fil_pos_pred, drop = FALSE]
    states_act_fil <- states_act_full[,fil_pos_act, drop = FALSE]

    states_pred_fil[states_pred_fil %in% c(0,1)] <- 0
    states_act_fil[states_act_fil %in% c(0,1)] <- 0

    ord <- cluster_split(states_pred_fil, states_act_fil, pos1 = pos1, pos2 = pos2)

    states_pred_ord <- states_pred_fil[ord$pred, , drop = FALSE]
    states_act_ord <- states_act_fil[ord$act, , drop = FALSE]

    pred_sm_df <- expand_grid(read = seq_len(nrow(states_pred_ord)),
                              pos = seq_len(ncol(states_pred_ord))) %>%
      mutate(state = factor(as.vector(t(states_pred_ord)), levels = c(-1, 0 ,2,3)), type = "pred")

    act_sm_df <- expand_grid(read = seq_len(nrow(states_act_ord)), pos = seq_len(ncol(states_act_ord))) %>%
      mutate(state = factor(as.vector(t(states_act_ord)), levels = c(-1, 0 ,2,3)), type = "act")

    sm_df <- rbind(pred_sm_df, act_sm_df)

    state_df <- bind_rows(tibble(
        state = names(ord$state_fraction_pred), frac = as.numeric(ord$state_fraction_pred), type = "pred"
      ), tibble(state = names(ord$state_fraction_act), frac = as.numeric(ord$state_fraction_act),
        type = "act")) %>%
      separate(state, into = c("left", "right"), sep = "_") %>%
      mutate(cat = case_when(left == "2" & right == "2" ~ "2-2",
        left == "2" & right != "2" ~ "2-non2", left != "2" & right == "2" ~ "non2-2", TRUE ~ "other"
      )) %>%
      group_by(type, cat) %>%
      summarise(frac = sum(frac), .groups = "drop")

    annot_df <- state_df %>%
      filter(cat != "other") %>%
      mutate(read_n = ifelse(type == "act", nrow(states_act_ord), nrow(states_pred_ord))) %>%
      group_by(type) %>%
      arrange(factor(cat, levels = c("2-2", "2-non2", "non2-2")), .by_group = TRUE) %>%
      mutate(ymin = lag(cumsum(frac), default = 0) * first(read_n),
        ymax = cumsum(frac) * first(read_n), ymid = (ymin + ymax) / 2, frac = frac*100) %>%
      ungroup()

    p_sm <- sm_df %>%
      ggplot() +
      facet_wrap(~type, scale = "free_y") +
      geom_raster(aes(x = pos, y = -read, fill = state)) +
      labs(x = "Genomic position", y = "Single molecules", fill = "State", title = paste("Pred= ", mi_pred, " Act=", mi_act)) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      scale_fill_manual(values = c("#FFFFFF", "#DEDEDE", "#149246", "#000000")) +
      geom_rect(data = annot_df, aes(ymin = -ymin, ymax = -ymax, xmin = 250, xmax = 270, color = cat)) +
      geom_text(data = annot_df, aes(label = round(frac), y = -ymid, x = 300), size = 5/2.8)

    pdf(file.path(output_dir, paste0(stem, "_single_molecules.pdf")), width = 20/2.54, height = 20/2.54)
    fix <- set_panel_size(p_sm, height = unit(1.5, "cm"), width = unit(2, "cm"))
    grid.arrange(fix)
    dev.off()

  }
}

main()
