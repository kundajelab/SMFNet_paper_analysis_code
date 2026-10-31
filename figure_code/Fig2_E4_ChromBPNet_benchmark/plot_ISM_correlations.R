#Run from the repository root: Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_ISM_correlations.R
#Run compute_ChromBPNet_ISM.py first.
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
output_dir <- file.path("output", "Fig2_E4_ChromBPNet_benchmark", "ISM_correlations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#theme is just for plotting
theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(),
                  plot.title = element_text(size = 6), axis.text = element_text(size = 6),
                  title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")))

ism_h5 <- file.path(input_dir, "SMFNet_endogenous_ISM.h5")
chrombpnet_h5 <- file.path("output", "Fig2_E4_ChromBPNet_benchmark", "ChromBPNet_endogenous_ISM.h5")
tfbs_file <- file.path(input_dir, "SMFNet_annotated_motif_hits_for_codependency_analysis.tsv")

#Process TFBS file
tfbs_df <- read.delim(tfbs_file) %>%
  filter(chr %in% c("chr1", "chr3", "chr6")) %>%
  select(chr, start, end, pattern, motif)

tfbs_gr <- GRanges(
  seqnames = tfbs_df$chr,
  ranges = IRanges(
    start = tfbs_df$start,
    end = tfbs_df$end
  ),
  pattern = tfbs_df$pattern,
  motif = tfbs_df$motif
)

smf_h5 <- H5File$new(ism_h5, mode = "r")
chrom_h5 <- H5File$new(chrombpnet_h5, mode = "r")

region_ids <- names(smf_h5[["scaled_ism_sum"]])

region_windows <- GRanges(
  seqnames = str_extract(region_ids, "^chr[^_]+"),
  ranges = IRanges(
    start = as.integer(str_match(region_ids, "_(\\d+)-")[, 2]) + 512,
    end = as.integer(str_match(region_ids, "-(\\d+)$")[, 2]) - 512
  )
)

region_ids <- region_ids[overlapsAny(region_windows, tfbs_gr)]

tfbs_list <- list()

for (i in seq_along(region_ids)) {

  rid <- region_ids[i]

  print(paste0("[", i, "/", length(region_ids), "] Processing ", rid, "..."))

  chr   <- str_extract(rid, "^chr[^_]+")
  window_start <- as.integer(str_match(rid, "_(\\d+)-")[,2]) + 512
  window_end   <- as.integer(str_match(rid, "-(\\d+)$")[,2]) - 512
  window_id <- paste0(chr, "_", window_start, "_", window_end)
  window_gr <- GRanges(seqnames = chr, ranges = IRanges(start = window_start, end = window_end))

  tf_hits <- findOverlaps(tfbs_gr, window_gr)

  scores <- tryCatch(
    smf_h5[["scaled_ism_sum"]][[rid]][],
    error = function(e) NULL
  )

  base_scores <- tryCatch(
    smf_h5[["ism_sum"]][[rid]][],
    error = function(e) NULL
  )

  chrom_scores <- tryCatch(
    chrom_h5[["ism_sum"]][[rid]][],
    error = function(e) NULL
  )

  tfbs_subset_df <- as.data.frame(tfbs_gr[queryHits(tf_hits)]) %>%
    mutate(center = (start + end) / 2, rel_start = start - window_start, rel_end = end - window_start) %>%
    filter(rel_start >= 1 & rel_end <= 1024) %>%
    unique()

  if (is.null(scores) || is.null(base_scores) || is.null(chrom_scores) || nrow(tfbs_subset_df) == 0) {
    message("Skipping missing region: ", rid)
    next
  }

  out_df <- tfbs_subset_df %>%
    rowwise() %>%
    mutate(smf_mean = mean(scores[rel_start:rel_end], na.rm = TRUE),
      base_smf_mean = mean(base_scores[rel_start:rel_end], na.rm = TRUE),
      chrom_mean = mean(chrom_scores[(rel_start+512):(rel_end+512)], na.rm = TRUE),
      smf_sum = sum(scores[rel_start:rel_end], na.rm = TRUE),
      base_smf_sum = sum(base_scores[rel_start:rel_end], na.rm = TRUE),
      chrom_sum = sum(chrom_scores[(rel_start+512):(rel_end+512)], na.rm = TRUE),
      chrom_max = which.max(chrom_scores[(rel_start+512):(rel_end+512)]),
      smf_max = which.max(scores[rel_start:rel_end]),
      base_smf_max = which.max(base_scores[rel_start:rel_end]),
      region_id = rid,
      window = window_id) %>%
    ungroup()

  tfbs_list[[rid]] = out_df

}

smf_h5$close_all()
chrom_h5$close_all()

attr_df <- do.call(rbind, tfbs_list)

total_cor <- cor(attr_df$base_smf_sum, attr_df$chrom_sum)

#Plot scatter
cor_attr <- attr_df %>%
  ggplot(aes(x = base_smf_sum, y = chrom_sum)) +
  geom_hex(bins = 200) +
  scale_fill_viridis_c(option = "A", trans = "log10") +
  annotate(geom = "text", x = -800, y = 15, size = 5 / 2.8, label = paste0("r=", round(total_cor, 2))) +
  theme(legend.position = "right", legend.text = element_text(size = 6)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed")

pdf(file.path(output_dir, "cor_ism_chrom_smf.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(cor_attr, height = unit(3, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
dev.off()

#Preserve the original motif panels: scaled scores on x, unscaled-score Pearson labels.
total_cor_motif <- attr_df %>%
  group_by(motif) %>%
  summarize(cor = cor(base_smf_sum, chrom_sum), min_smf = min(base_smf_sum), max_chrom = max(chrom_sum)) %>%
  mutate(y = max_chrom - 0.1 * max_chrom,x = min_smf - 0.2 * min_smf, label = paste0("r = ", round(cor, 2)))

cor_attr_motif <- attr_df %>%
  ggplot(aes(x = smf_sum, y = chrom_sum)) +
  facet_wrap(~motif, scales = "free") +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(option = "A", trans = "log10") +
  geom_text(data = total_cor_motif, inherit.aes = FALSE,
            aes(label = label, x = x, y = y), size = 5 / 2.8) +
  theme(legend.position = "right", legend.text = element_text(size = 6)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = 0), linetype = "dashed")

pdf(file.path(output_dir, "cor_ism_chrom_smf_motifs.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(cor_attr_motif, height = unit(1, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
dev.off()
