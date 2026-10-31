#Run from the repository root: Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_marginalizations.R
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(egg)
library(gridExtra)
library(hdf5r)

input_dir <- file.path("input_files", "Fig2_E4_ChromBPNet_benchmark")
output_dir <- file.path("output", "Fig2_E4_ChromBPNet_benchmark", "marginalizations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(),
                  plot.title = element_text(size = 6), axis.text = element_text(size = 6),
                  title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")))

#Load ChromBPNet profiles for the total marginalization correlation.
chrom_df <- read.delim(file.path(input_dir, "ChromBPNet_marginalization_profiles.tsv")) %>%
  separate(region_id, c("region_id", "motif"), sep = "__") %>%
  mutate(pos = pos + 1, model = "ChromBPNet")

#Compare central accessibility with the original ChromBPNet profile sum.
#Keep one row per selected position: this preserves the original weighted motif median.
smf_states <- file.path(input_dir, "SMFNet_marginalization_states.h5")
states_h5 <- H5File$new(smf_states, mode = "r")

#State-call files retain their original /smf_states and /smf_pos groups.
region_list <- states_h5[["smf_states"]]$ls()$name

extract_states <- function(states_h5, region, center_pos = 511.5, tfbs_dist = 25) {

  smf_pos <- states_h5[["smf_pos"]][[region]][]
  states_mat <- t(states_h5[["smf_states"]][[region]]$read(drop = FALSE))

  dists <- abs(smf_pos - center_pos)

  if (min(dists) > tfbs_dist) {
    return(NULL)
  }

  keep_idx <- dists < tfbs_dist

  states_vec <- states_mat[, keep_idx, drop = FALSE]

  oc <- sum(apply(states_vec, 1, \(x) any(x == 3))) / nrow(states_vec)
  ac <- 1 - oc

  data.frame(region = region,
             ac_perc = ac,
             oc_perc = oc,
             tfbs_pos = smf_pos[keep_idx])
}

smf_tfbs_list <- vector("list", length(region_list))

for (i in seq_along(region_list)) {
  r <- region_list[i]
  if (i == 1L || i %% 250L == 0L) message(i, "/", length(region_list), ": ", r)
  smf_tfbs_list[[i]] <- extract_states(states_h5, r)
}

states_h5$close_all()

tfbs_df <- do.call(rbind, smf_tfbs_list)
smf_perc <- tfbs_df %>%
  separate(region, c("region_id", "motif"), sep = "__") %>%
  group_by(motif) %>%
  summarize(median_ac = median(ac_perc))

chrom_sum <- chrom_df %>%
  group_by(motif, pos, model) %>%
  summarize(signal = median(signal)) %>%
  filter(between(pos, 392, 542)) %>%
  ungroup() %>%
  group_by(motif) %>%
  summarize(med_chrom = sum(signal))

combi_df <- left_join(smf_perc, chrom_sum) %>%
  filter(!str_detect(motif, "alt|double|pattern_19_unknown"))

cor_margs <- cor(combi_df$median_ac, combi_df$med_chrom)

plot_cor <- combi_df %>%
  ggplot(aes(y = med_chrom, x = median_ac)) +
  geom_point() +
  geom_point(data = combi_df[combi_df$motif == "noMotif",], shape = 21, fill = "#676767") +
  geom_text(aes(label = motif), size = 4 / 2.8) +
  annotate(geom = "text", label = paste0("r=", round(cor_margs, 2)), size = 5 / 2.8, x = 0.3, y = 225) +
  scale_y_continuous(limits = c(20, 250)) +
  scale_x_continuous(limits = c(0.2, 0.9))

pdf(file.path(output_dir, "full_chrom_smfnet_comp.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(plot_cor, height = unit(3, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
dev.off()
