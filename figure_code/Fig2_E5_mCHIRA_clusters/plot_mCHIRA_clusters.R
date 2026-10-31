#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig2_E5_mCHIRA_clusters/plot_mCHIRA_clusters.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(egg)
library(gridExtra)
library(hdf5r)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(),
                  plot.title = element_text(size = 6), axis.text = element_text(size = 6),
                  title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), plot.margin = grid::unit(c(0, 0, 0, 0), "cm")))

PLOT_DIR <- "output/Fig2_E5_mCHIRA_clusters/"

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
cluster_h5 <- "output/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_clusters.h5"
ACT_H5 <- "input_files/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_Experiment_reads.h5"
PRED_H5 <- "input_files/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_SMFNet_reads.h5"

#Compare cluster fractions and changes relative to wildtype.
df <- read.delim("output/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_cluster_fractions.tsv")

tfbs_info <- data.frame(motif = c("YY1", "NRF1", "KLF", "KLF", "MAX-MYC", "KLF"),
                        center = c(64.5, 100.0, 143.5, 155.5, 173.5, 196.0))

new_cluster <- c("2" = "A", "1" = "B", "3" = "C", "0" = "D", "4" = "E", "5" = "F", "6" = "G", "7" = "H")

fil_df <- df %>%
  dplyr::filter(grepl("6005", group))

fil_df$cluster <- new_cluster[as.character(fil_df$cluster)]

wide_df <- fil_df %>%
  dplyr::select(-count) %>%
  tidyr::pivot_wider(names_from = sample, values_from = fraction)

cluster_cor <- cor(wide_df$act, wide_df$pred)

cluster_scatter <- wide_df %>%
  ggplot() +
  geom_point(aes(x = act, y = pred, color = factor(cluster)), alpha = 0.5, size = 0.8) +
  annotate(geom = "text", label = paste0("r=", round(cluster_cor, 2)), x = 0.1, y = 0.55, size = 5/2.8) +
  coord_cartesian(xlim = c(0, 0.65), ylim = c(0, 0.65)) +
  scale_color_manual(values = c("#4c134c", "#7f277b", "#88439a", "#8c6caf", "#8c96c6", "#9fbcda", "#bfd4e6", "#e0ebf3")) +
  geom_smooth(aes(x = act, y = pred), color = "black", se = FALSE, method = "lm", size = 0.5)

pdf(paste0(PLOT_DIR, "rmce_cluster_scatter.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(cluster_scatter, height = grid::unit(3, "cm"), width = grid::unit(3, "cm"))
grid.arrange(fix)
dev.off()

wt_bar <- fil_df %>%
  dplyr::filter(group == "CRE_6005") %>%
  ggplot() +
  geom_bar(aes(y = cluster, x = fraction, fill = sample), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#ABABAB","#636ab1")) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(0, 0.6))

pdf(paste0(PLOT_DIR, "rmce_wt_bar_freq.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(wt_bar, height = grid::unit(2.5, "cm"), width = grid::unit(1.5, "cm"))
grid.arrange(fix)
dev.off()

ref <- fil_df[fil_df$group == "CRE_6005", c("sample", "cluster", "fraction")]
names(ref)[3] <- "ref_fraction"

delta_df <- merge(
  fil_df[fil_df$group != "CRE_6005", ],
  ref,
  by = c("sample", "cluster")
)

delta_df$delta <- delta_df$fraction - delta_df$ref_fraction

plot_delta_df <- delta_df %>%
  dplyr::select(-c(count, ref_fraction, fraction)) %>%
  tidyr::pivot_wider(names_from = sample, values_from = delta)

delta_cor <- cor(plot_delta_df$act, plot_delta_df$pred)

delta_plot <- plot_delta_df %>%
  ggplot() +
  geom_point(aes(x = act, y = pred, color = cluster),  alpha = 0.5, size = 0.8) +
  coord_cartesian(xlim = c(-0.45, 0.25), ylim = c(-0.45, 0.25)) +
  scale_color_manual(values = c("#4c134c", "#7f277b", "#88439a", "#8c6caf", "#8c96c6", "#9fbcda", "#bfd4e6", "#e0ebf3")) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  annotate(geom = "text", label = paste0("r=", round(delta_cor, 2)), x = -0.4, y = 0.2, size = 5/2.8)

pdf(paste0(PLOT_DIR, "rmce_delta_scatter.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(delta_plot, height = grid::unit(3, "cm"), width = grid::unit(3, "cm"))
grid.arrange(fix)
dev.off()

out_file <- wide_df %>%
  dplyr::rename(act_frac = "act", pred_frac = "pred") %>%
  left_join(plot_delta_df) %>%
  dplyr::rename(act_delta = "act", pred_delta = "pred")

write.table(out_file, "output/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_cluster_perturbations.tsv",
            sep = "\t", quote = FALSE)

delta_nrf1_bar <- delta_df %>%
  dplyr::filter(group == "CRE_6005_mutant_1") %>%
  ggplot() +
  facet_wrap(~sample) +
  geom_bar(aes(y = cluster, x = delta, fill = cluster), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#4c134c", "#7f277b", "#88439a", "#8c6caf", "#8c96c6", "#9fbcda", "#bfd4e6", "#e0ebf3")) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(-0.4, 0.2)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed")

pdf(paste0(PLOT_DIR, "rmce_dnrf1_delta.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(delta_nrf1_bar, height = grid::unit(2, "cm"), width = grid::unit(1.75, "cm"))
grid.arrange(fix)
dev.off()

delta_myc_bar <- delta_df %>%
  dplyr::filter(group == "CRE_6005_mutant_4") %>%
  ggplot() +
  facet_wrap(~sample) +
  geom_bar(aes(y = cluster, x = delta, fill = cluster), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#4c134c", "#7f277b", "#88439a", "#8c6caf", "#8c96c6", "#9fbcda", "#bfd4e6", "#e0ebf3")) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(-0.4, 0.2)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed")

pdf(paste0(PLOT_DIR, "rmce_dmyc_delta.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(delta_myc_bar, height = grid::unit(2, "cm"), width = grid::unit(1.75, "cm"))
grid.arrange(fix)
dev.off()

h5_cl   <- H5File$new(cluster_h5, mode = "r")
h5_act  <- H5File$new(ACT_H5, mode = "r")
h5_pred <- H5File$new(PRED_H5, mode = "r")

orig_h5 <- list(
  act = h5_act,
  pred = h5_pred
)

region_paths <- lapply(orig_h5, function(h5) {
  setNames(sprintf("regions/%08d", seq_along(h5[["region_ids"]][]) - 1L),
           h5[["region_ids"]][])
})
roles <- c(act = "observed", pred = "predicted")

read_list <- list()

for (sample in c("act", "pred")) {

  groups <- names(h5_cl[[sample]])

  groups <- groups[grepl("^CRE_6005($|_mutant_)", groups)]

  for (group in groups) {

    cl <- h5_cl[[sample]][[group]][["cluster_labels"]][]
    read_ids <- h5_cl[[sample]][[group]][["read_ids"]][]

    stopifnot(length(cl) == length(read_ids))

    read_group <- orig_h5[[sample]][[region_paths[[sample]][[group]]]][[roles[[sample]]]]
    mat <- read_group[["smf"]][,]
    pos <- read_group[["positions"]][]

    if (nrow(mat) != length(pos)) {
      mat <- t(mat)
    }

    idx <- read_ids + 1L

    if (max(idx) > ncol(mat)) {
      stop(
        sample, " / ", group,
        ": cluster read ID exceeds number of original reads"
      )
    }

    mat_sel <- mat[, idx, drop = FALSE]

    tmp <- data.frame(
      pos = rep(pos, times = ncol(mat_sel)),
      state = as.vector(mat_sel),
      read_index = rep(seq_len(ncol(mat_sel)), each = length(pos)),
      cluster = rep(cl, each = length(pos))
    )

    tmp$sample <- sample
    tmp$group <- group

    tmp$molecule <- paste(
      sample,
      group,
      read_ids[tmp$read_index],
      sep = "_"
    )

    read_list[[paste(sample, group, sep = "_")]] <- tmp
  }
}

h5_cl$close_all()
h5_act$close_all()
h5_pred$close_all()

original_cluster_reads <- dplyr::bind_rows(read_list)

original_cluster_reads$cluster <- new_cluster[
  as.character(original_cluster_reads$cluster)
]

original_cluster_reads$cluster <- factor(
  original_cluster_reads$cluster,
  levels = LETTERS[1:8]
)

wt_df_all <- original_cluster_reads %>%
  filter(group == "CRE_6005") %>%
  filter(state != -1) %>%
  group_by(pos, sample) %>%
  summarize(signal = mean(state))

wt_avg <- wt_df_all %>%
  ggplot(aes(x = pos, y = signal, color = sample)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#ABABAB","#636ab1")) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 251)) +
  geom_rect(data = tfbs_info,
            aes(xmin = center - 5, xmax = center + 5,ymin = 0, ymax = 0.1), inherit.aes = FALSE, fill = "#DEDEDE") +
  geom_text(data = tfbs_info,
            aes(x = center, y = 0.1, label = motif), inherit.aes = FALSE, size = 5/2.8) +
  labs(x = "Position", y = "Mean SMF signal (%)", color = NULL)

pdf(paste0(PLOT_DIR, "rmce_wt_avg.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(wt_avg, height = grid::unit(2, "cm"), width = grid::unit(4, "cm"))
grid.arrange(fix)
dev.off()

cluster_sum <- original_cluster_reads %>%
  filter(state != -1) %>%
  group_by(pos, sample, cluster) %>%
  summarize(signal = mean(state))

cluster_smf_signal <- cluster_sum %>%
  ggplot(aes(x = pos, y = signal, color = sample)) +
  facet_wrap(~cluster) +
  geom_line() +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("#ABABAB","#636ab1")) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_rect(data = tfbs_info,
            aes(xmin = center - 5, xmax = center + 5,ymin = 0, ymax = 0.1), inherit.aes = FALSE, fill = "#DEDEDE") +
  geom_text(data = tfbs_info,
            aes(x = center, y = 0.1, label = motif), inherit.aes = FALSE, size = 5/2.8) +
  labs(x = "Position", y = "Mean SMF signal (%)", color = NULL)

pdf(paste0(PLOT_DIR, "rmce_all_cluster_smf_avg.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(cluster_smf_signal, height = grid::unit(1.5, "cm"), width = grid::unit(2, "cm"))
grid.arrange(fix)
dev.off()

cluster_sum_wt <- original_cluster_reads %>%
  filter(state != -1) %>%
  filter(group == "CRE_6005") %>%
  group_by(pos, sample, cluster) %>%
  summarize(signal = mean(state))

cluster_smf_signal_wt <- cluster_sum_wt %>%
  ggplot(aes(x = pos, y = signal, color = sample)) +
  facet_wrap(~cluster) +
  geom_line() +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("#ABABAB","#636ab1")) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_rect(data = tfbs_info,
            aes(xmin = center - 5, xmax = center + 5,ymin = 0, ymax = 0.1), inherit.aes = FALSE, fill = "#DEDEDE") +
  geom_text(data = tfbs_info,
            aes(x = center, y = 0.1, label = motif), inherit.aes = FALSE, size = 5/2.8) +
  labs(x = "Position", y = "Mean SMF signal (%)", color = NULL)

pdf(paste0(PLOT_DIR, "rmce_wt_cluster_smf_avg.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(cluster_smf_signal_wt, height = grid::unit(1.5, "cm"), width = grid::unit(4, "cm"))
grid.arrange(fix)
dev.off()

h5_cl <- H5File$new(cluster_h5, mode = "r")

cluster_read_list <- list()

for (sample in c("act", "pred")) {

  groups <- names(h5_cl[[sample]])

  groups <- groups[grepl("^CRE_6005($|_mutant_)", groups)]

  for (group in groups) {

    mat <- h5_cl[[sample]][[group]][["expanded_reads"]][,]
    cl <- h5_cl[[sample]][[group]][["cluster_labels"]][]
    read_ids <- h5_cl[[sample]][[group]][["read_ids"]][]

    stopifnot(length(cl) == length(read_ids))

    if (nrow(mat) != length(cl)) {
      mat <- t(mat)
    }

    stopifnot(nrow(mat) == length(cl))

    keep <- !is.na(cl) & cl >= 0

    mat <- mat[keep, , drop = FALSE]
    cl <- cl[keep]
    read_ids <- read_ids[keep]

    tmp <- data.frame(
      sample = sample,
      group = group,
      read_id = rep(read_ids, each = ncol(mat)),
      molecule = rep(
        paste(sample, group, read_ids, sep = "_"),
        each = ncol(mat)
      ),
      cluster = rep(cl, each = ncol(mat)),
      pos = rep(0:(ncol(mat) - 1), times = nrow(mat)),
      state = as.vector(t(mat))
    )

    cluster_read_list[[paste(sample, group, sep = "_")]] <- tmp
  }
}

h5_cl$close_all()

cluster_single_reads <- dplyr::bind_rows(cluster_read_list)

cluster_single_reads$cluster <- new_cluster[
  as.character(cluster_single_reads$cluster)
]

cluster_single_reads$cluster <- factor(
  cluster_single_reads$cluster,
  levels = LETTERS[1:8]
)

read_order <- cluster_single_reads %>%
  group_by(cluster, read_id, group, molecule, sample) %>%
  summarize(avg = mean(state), .groups = "drop") %>%
  group_by(cluster, group, sample) %>%
  arrange(avg, .by_group = TRUE) %>%
  mutate(read_order = row_number())

sm_plot_act <- cluster_single_reads %>%
  left_join(read_order) %>%
  filter(sample == "act") %>%
  filter(group == "CRE_6005") %>%
  ggplot(aes(x = as.numeric(pos), y = -read_order, fill = as.factor(state))) +
  facet_wrap(~cluster) +
  geom_raster() +
  scale_y_discrete(name = "Single molecules (Nucleosomes)", expand = c(0, 0)) +
  scale_fill_manual(values = c("#DEDEDE", "#000000"), name = "Nuc State") +
  theme_void() +
  theme(axis.text = element_blank(), legend.position = "none",
        axis.title.y = element_blank(), text = element_text(size = 6)) +
  scale_x_continuous(name = "Relative position", expand = c(0, 0))

pdf(paste0(PLOT_DIR, "wt_cluster_sm_act.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(sm_plot_act, height = grid::unit(1, "cm"), width = grid::unit(3, "cm"))
grid.arrange(fix)
dev.off()

sm_plot_pred <- cluster_single_reads %>%
  left_join(read_order) %>%
  filter(sample == "pred") %>%
  filter(group == "CRE_6005") %>%
  ggplot(aes(x = as.numeric(pos), y = -read_order, fill = as.factor(state))) +
  facet_wrap(~cluster) +
  geom_raster() +
  scale_y_discrete(name = "Single molecules (Nucleosomes)", expand = c(0, 0)) +
  scale_fill_manual(values = c("#DEDEDE", "#000000"), name = "Nuc State") +
  theme_void() +
  theme(axis.text = element_blank(), legend.position = "none",
        axis.title.y = element_blank(), text = element_text(size = 6)) +
  scale_x_continuous(name = "Relative position", expand = c(0, 0))

pdf(paste0(PLOT_DIR, "wt_cluster_sm_pred.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(sm_plot_pred, height = grid::unit(1, "cm"), width = grid::unit(3, "cm"))
grid.arrange(fix)
dev.off()

read_order_all <- cluster_single_reads %>%
  group_by(cluster, read_id,  molecule, sample) %>%
  summarize(avg = mean(state), .groups = "drop") %>%
  group_by(cluster, sample) %>%
  arrange(avg, .by_group = TRUE) %>%
  mutate(read_order = row_number())

sm_plot_act_all <- cluster_single_reads %>%
  left_join(read_order_all) %>%
  filter(sample == "act") %>%
  ggplot(aes(x = as.numeric(pos), y = -read_order, fill = as.factor(state))) +
  facet_wrap(~cluster) +
  geom_raster() +
  scale_y_discrete(name = "Single molecules (Nucleosomes)", expand = c(0, 0)) +
  scale_fill_manual(values = c("#DEDEDE", "#000000"), name = "Nuc State") +
  theme_void() +
  theme(axis.text = element_blank(), legend.position = "none",
        axis.title.y = element_blank(), text = element_text(size = 6)) +
  scale_x_continuous(name = "Relative position", expand = c(0, 0))

pdf(paste0(PLOT_DIR, "all_cluster_sm_act.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(sm_plot_act_all, height = grid::unit(1, "cm"), width = grid::unit(3, "cm"))
grid.arrange(fix)
dev.off()

sm_plot_pred_all <- cluster_single_reads %>%
  left_join(read_order_all) %>%
  filter(sample == "pred") %>%
  ggplot(aes(x = as.numeric(pos), y = -read_order, fill = as.factor(state))) +
  facet_wrap(~cluster) +
  geom_raster() +
  scale_y_discrete(name = "Single molecules (Nucleosomes)", expand = c(0, 0)) +
  scale_fill_manual(values = c("#DEDEDE", "#000000"), name = "Nuc State") +
  theme_void() +
  theme(axis.text = element_blank(), legend.position = "none",
        axis.title.y = element_blank(), text = element_text(size = 6)) +
  scale_x_continuous(name = "Relative position", expand = c(0, 0))

pdf(paste0(PLOT_DIR, "all_cluster_sm_pred.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(sm_plot_pred_all, height = grid::unit(1, "cm"), width = grid::unit(3, "cm"))
grid.arrange(fix)
dev.off()
