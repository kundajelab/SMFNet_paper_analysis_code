#Run from the repository root: Rscript figure_code/Fig1_E2_E3_Extra_models_benchmarking/plot_Sonmezer_qc.R
library(ggplot2)
library(dplyr)
library(tidyr)
library(egg)
library(gridExtra)
library(hdf5r)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6),
                  plot.title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  title = element_text(size = 6), legend.position = "none",
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")))

input_dir <- file.path("input_files", "Fig1_E2_E3_Extra_models_benchmarking")
output_dir <- file.path("output", "Fig1_E2_E3_Extra_models_benchmarking")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
#Read smfnet matrices as molecules by positions.
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

pred_avg_sonmezer <- file.path(input_dir, "avg_TestRegions_Sonmezer_SMFNet.tsv")
pred_avg_ont <- file.path(input_dir, "avg_TestRegions_SMFNet.tsv")
pred_handle <- H5File$new(file.path(input_dir, "mESC_Sonmezer_testRegions_SMFNet_reads.h5"), mode = "r")
act_handle <- H5File$new(file.path(input_dir, "mESC_Sonmezer_testRegions_Experiment_reads.h5"), mode = "r")
pred_index <- sample_index(pred_handle)
act_index <- sample_index(act_handle)

#Expand comma-separated site profiles, then match models by region and position.
avg_sonmezer_df <- read.delim(pred_avg_sonmezer,
                      header = TRUE) %>%
  separate_longer_delim(c(smf_pos, avg_signal), delim = ",") %>%
  mutate(type = "sonmezer")

avg_ont_df <- read.delim(pred_avg_ont,
                              header = TRUE) %>%
  separate_longer_delim(c(smf_pos, avg_signal), delim = ",") %>%
  mutate(type = "ont")

pred_avg_total <- rbind(avg_sonmezer_df, avg_ont_df) %>%
  pivot_wider(names_from = type, values_from = avg_signal) %>%
  na.omit()

#Compare predicted averages at matching sites
total_cor <- cor(as.numeric(pred_avg_total$sonmezer), as.numeric(pred_avg_total$ont))

scatter_cor <- pred_avg_total %>%
  ggplot(aes(y = as.numeric(sonmezer), x = as.numeric(ont))) +
  geom_hex(bins = 500) +
  scale_fill_viridis_c(option = "A", trans = "log10") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  annotate(geom = "text", x = 0.15, y = 0.95, size = 5 / 2.8, label = paste0("r=", round(total_cor, 2))) +
  theme(legend.position = "right", legend.text = element_text(size = 6))

pdf(file.path(output_dir, "TestRegions_smfpos_scatter_sonmezer_ont.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(scatter_cor, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
dev.off()

#Use regions with at least 70 shared sites for the observed-data comparisons.
loop_regions <- pred_avg_total %>%
  group_by(region) %>%
  tally() %>%
  filter(n >= 70) %>%
  pull(region) %>%
  unique()

#Some test regions have no observed short-read data.
loop_regions <- intersect(loop_regions, names(act_index))

#Align observed and predicted sites in the central 1024 bp window.
subset_to_cols <- function(x, cols) {
  mat <- x$mat[, match(cols, x$pos), drop = FALSE]
  colnames(mat) <- cols
  mat
}

#Compare predicted and observed co-variation using pairwise log2 odds ratios.
#A 0.5 correction keeps odds ratios finite when a contingency cell is empty.

pairwise_odds_ratios <- function(act_mat, pred_mat, min_pair_cov = 50,
                                 invalid_id = -1, correction = 0.5) {

  pair_counts <- function(mat) {
    observed <- mat != invalid_id
    is_zero  <- mat == 0
    is_one   <- mat == 1

    list(
      coverage = crossprod(observed),
      n00 = crossprod(is_zero),
      n01 = crossprod(is_zero, is_one),
      n10 = crossprod(is_one, is_zero),
      n11 = crossprod(is_one)
    )
  }

  act  <- pair_counts(act_mat)
  pred <- pair_counts(pred_mat)

  #Same ordering as combn(seq_len(ncol(mat)), 2)
  pair_idx <- which(upper.tri(act$coverage), arr.ind = TRUE)

  act_cov  <- act$coverage[pair_idx]
  pred_cov <- pred$coverage[pair_idx]

  keep <- act_cov >= min_pair_cov & pred_cov > 0

  if (!any(keep)) {
    return(data.frame())
  }

  pair_idx <- pair_idx[keep, , drop = FALSE]
  act_cov  <- act_cov[keep]
  pred_cov <- pred_cov[keep]

  get_counts <- function(x) {
    list(
      n00 = x$n00[pair_idx] + correction,
      n01 = x$n01[pair_idx] + correction,
      n10 = x$n10[pair_idx] + correction,
      n11 = x$n11[pair_idx] + correction
    )
  }

  act_count  <- get_counts(act)
  pred_count <- get_counts(pred)

  act_log2_or <- log2(
    (act_count$n00 * act_count$n11) /
      (act_count$n01 * act_count$n10)
  )

  pred_log2_or <- log2(
    (pred_count$n00 * pred_count$n11) /
      (pred_count$n01 * pred_count$n10)
  )

  data.frame(
    pos1 = as.numeric(colnames(act_mat)[pair_idx[, 1]]),
    pos2 = as.numeric(colnames(act_mat)[pair_idx[, 2]]),
    act_cov = act_cov,
    pred_cov = pred_cov,
    act_log2_or = act_log2_or,
    pred_log2_or = pred_log2_or,
    log2_or_diff = pred_log2_or - act_log2_or
  )
}

odds_list <- list()
avg_list <- list()
min_reads <- 100
min_cols <- 15

for (i in seq_along(loop_regions)) {

  r <- loop_regions[i]

  print(paste0(i, "/", length(loop_regions), " processing ", r))
  start_loop <- Sys.time()

  sonmezer <- read_sample_region(pred_handle, pred_index, r, "predicted")
  act <- read_sample_region(act_handle, act_index, r, "observed")
  #Observed positions span 2048 bp; predictions refer to the central 1024 bp.
  act$pos <- act$pos - 512

  common_cols <- intersect(sonmezer$pos, act$pos)

  pred_sonmezer_mat <- subset_to_cols(sonmezer, common_cols)
  act_sonmezer_mat <- subset_to_cols(act, common_cols)

  if (nrow(act_sonmezer_mat) < min_reads) {
    cat(sprintf("   skipping: only %d act reads\n\n", nrow(act_sonmezer_mat)))
    next
  }

  if (ncol(act_sonmezer_mat) < min_cols) {
    cat(sprintf("   skipping: only %d common positions\n\n", ncol(act_sonmezer_mat)))
    next
  }

  #calc paired odds
  region_odds <- pairwise_odds_ratios(
    act_mat = act_sonmezer_mat,
    pred_mat = pred_sonmezer_mat,
    min_pair_cov = min_reads
  )

  if (nrow(region_odds) == 0) {
    cat("   skipping: no position pairs with sufficient coverage\n\n")
    next
  }

  region_odds$region <- r

  odds_list[[length(odds_list) + 1]] <- region_odds

  #Calculate per-position averages
  act_cov  <- colSums(act_sonmezer_mat != -1)
  pred_cov <- colSums(pred_sonmezer_mat != -1)

  #Average-signal plots require at least 100 calls at each site in both datasets.
  keep_pos <- act_cov >= min_reads & pred_cov >= min_reads

  if (any(keep_pos)) {

    act_avg  <- colSums(act_sonmezer_mat == 1) / act_cov
    pred_avg <- colSums(pred_sonmezer_mat == 1) / pred_cov

    avg_list[[length(avg_list) + 1]] <- data.frame(
      region = r,
      smf_pos = as.numeric(common_cols[keep_pos]),
      act_cov = act_cov[keep_pos],
      pred_cov = pred_cov[keep_pos],
      act_avg = act_avg[keep_pos],
      pred_avg = pred_avg[keep_pos],
      avg_diff = pred_avg[keep_pos] - act_avg[keep_pos]
    )
  }

  cat(sprintf(
    "   loop time: %.2f sec\n\n",
    as.numeric(difftime(Sys.time(), start_loop, units = "secs"))
  ))
}

odds_long <- do.call(rbind, odds_list)
avg_long <- do.call(rbind, avg_list)

write.table(odds_long, file.path(output_dir, "pos_odds_preds_vs_act.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(avg_long, file.path(output_dir, "pos_avg_preds_vs_act.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#Plot co-binding comp vs act
cor_odds <- cor(odds_long$act_log2_or, odds_long$pred_log2_or)

scatter_cor_odds <- odds_long %>%
  ggplot(aes(y = as.numeric(act_log2_or), x = as.numeric(pred_log2_or))) +
  geom_hex(bins = 500) +
  scale_fill_viridis_c(option = "A", trans = "log10") +
  annotate(geom = "text", x = -5, y = 12.5, size = 5 / 2.8, label = paste0("r=", round(cor_odds, 2))) +
  theme(legend.position = "right", legend.text = element_text(size = 6)) +
  coord_cartesian(xlim = c(-7.5, 15), ylim = c(-7.5, 15)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed")

pdf(file.path(output_dir, "TestRegions_odds_scatter_Sonmezer_predicted_vs_observed.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(scatter_cor_odds, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
dev.off()

#Plot co-binding by distance
odds_dist <- odds_long %>%
  mutate(dist = abs(pos1-pos2)) %>%
  pivot_longer(c(act_log2_or, pred_log2_or), names_to = "type", values_to = "or")

#Overlay the median odds ratio in each 1 bp distance bin.
smooth_df <- odds_dist %>%
  filter(dist >= 0, dist <= 300, or >= -7.5, or <= 15) %>%
  mutate(dist_bin = floor(dist)) %>%
  group_by(type, dist_bin) %>%
  summarize(dist = mean(dist), or = median(or), .groups = "drop")

dist_odds <- odds_dist %>%
  ggplot(aes(y = or, x = dist)) +
  facet_wrap(~type, ncol = 1) +
  geom_hex(bins = 500) +
  scale_fill_viridis_c(option = "A", trans = "log10") +
  geom_line(data = smooth_df, color = "black", linewidth = 0.5) +
  theme(legend.position = "right", legend.text = element_text(size = 6)) +
  coord_cartesian(xlim = c(0, 300), ylim = c(-7.5, 15)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed")

pdf(file.path(output_dir, "TestRegions_smfpos_dist_odds_sonmezer_pred_vs_act.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(dist_odds, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
dev.off()

#Calculate average correlation

total_cor <- cor(as.numeric(avg_long$act_avg), as.numeric(avg_long$pred_avg))

scatter_cor_act <- avg_long %>%
  ggplot(aes(y = as.numeric(act_avg), x = as.numeric(pred_avg))) +
  geom_hex(bins = 400) +
  scale_fill_viridis_c(option = "A", trans = "log10") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  annotate(geom = "text", x = 0.15, y = 0.95, size = 5 / 2.8, label = paste0("r=", round(total_cor, 2))) +
  theme(legend.position = "right", legend.text = element_text(size = 6))

pdf(file.path(output_dir, "TestRegions_smfpos_scatter_sonmezer_pred_act.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(scatter_cor_act, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
dev.off()

pred_handle$close_all()
act_handle$close_all()

