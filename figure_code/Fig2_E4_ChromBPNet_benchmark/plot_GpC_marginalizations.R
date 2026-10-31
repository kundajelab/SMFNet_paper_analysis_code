#Run from the repository root: Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_GpC_marginalizations.R
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(egg)
library(gridExtra)
library(hdf5r)
library(zoo)

input_dir <- file.path("input_files", "Fig2_E4_ChromBPNet_benchmark")
output_dir <- file.path("output", "Fig2_E4_ChromBPNet_benchmark", "GpC_marginalizations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#Eight precomputed marginalization runs: JASPAR motifs and RFX3 for each cell type.
RUNS <- c("jaspar_mESC", "jaspar_NPC", "jaspar_C2C12", "jaspar_MEL",
          "rfx3_mESC", "rfx3_NPC", "rfx3_C2C12", "rfx3_MEL")
FILES <- setNames(file.path(input_dir, paste0("GpC_", sub("^(.*)_(.*)$", "\\2_\\1", RUNS),
                                            "_marginalizations.h5")), RUNS)

theme_set(theme_classic() +
            theme(
              text = element_text(size = 6),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
              axis.line = element_blank(),
              strip.background = element_blank(),
              legend.position = "bottom"
            ))

#Average reads within each background, then take the median across backgrounds.
median_profile <- function(dataset, n_backgrounds, motif_index = NULL) {
  dims <- dataset$dims
  profiles <- lapply(seq_len(n_backgrounds), function(background) {
    #hdf5r reverses stored dimensions: positions, reads, backgrounds, motifs.
    indices <- list(seq_len(dims[1]), seq_len(dims[2]), background)
    if (!is.null(motif_index)) indices <- c(indices, list(motif_index))
    values <- dataset$read(args = indices, drop = FALSE)
    mat <- matrix(values, nrow = dims[1], ncol = dims[2])
    mat[mat == -1] <- NA_real_
    rowMeans(mat, na.rm = TRUE)
  })
  apply(do.call(cbind, profiles), 1, median, na.rm = TRUE)
}

read_profiles <- function(path, cell_type) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  stopifnot(h5attr(h5, "schema") == "smfnet.marginalizations",
            h5attr(h5, "schema_version") == 1)
  motifs <- as.character(h5[["motif_names"]][])
  n_backgrounds <- as.integer(h5attr(h5, "num_backgrounds"))
  background <- median_profile(h5[["samples/background"]], n_backgrounds)

  bind_rows(lapply(seq_along(motifs), function(m) {
    message(cell_type, ": ", motifs[m])
    profiles <- list(
      Background = background,
      Control = median_profile(h5[["samples/control"]], n_backgrounds, m),
      Motif = median_profile(h5[["samples/motif"]], n_backgrounds, m)
    )
    bind_rows(lapply(names(profiles), function(condition) {
      signal <- profiles[[condition]]
      tibble(cell_type = cell_type, motif = motifs[m], condition = condition,
             pos = seq_along(signal) - (length(signal) + 1) / 2, signal = signal)
    }))
  }))
}

profiles <- bind_rows(lapply(names(FILES), function(cell_type) {
  read_profiles(FILES[[cell_type]], cell_type)
})) %>%
  group_by(cell_type, motif, condition) %>%
  arrange(pos, .by_group = TRUE) %>%
  mutate(signal = zoo::rollapply(
    signal, width = 10, FUN = function(x) mean(x, na.rm = TRUE),
    align = "center", fill = NA, partial = TRUE
  )) %>%
  ungroup()

plot_motifs <- c( "NRF1(+)", "OCT4::SOX2(+)", "MAF::NFE2(+)", "GATA1::TAL1(+)", "RFX3(+)")

marg_plot <- profiles %>%
  filter(between(pos, -200, 200), str_ends(motif, fixed("(+)"))) %>%
  filter(condition != "Control") %>%
  filter(motif %in% plot_motifs) %>%
  mutate(cell_type = str_remove(cell_type, "^(jaspar|rfx3)_")) %>%
  mutate(cell_type = factor(cell_type, levels = c("mESC", "C2C12",  "MEL", "NPC")),
         motif = factor(motif, levels = plot_motifs)) %>%
  ggplot(aes(pos, signal, color = condition)) +
  geom_line(linewidth = 0.5) +
  facet_grid(cell_type ~ motif) +
  scale_color_manual(values = c(Background = "#ABABAB", Motif = "#636ab1")) +
  scale_x_continuous(limits = c(-200, 200), breaks = c(-150, 0, 150)) +
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7)) +
  labs(x = "Relative position (bp)",
       y = "GpC-only SMFNet signal (median, %)", color = NULL)

pdf(file.path(output_dir, "celltype_gc_marginalizations.pdf"),
      width = 20/2.54, height = 20/2.54)
  grid.arrange(set_panel_size(marg_plot, height = unit(1.25, "cm"), width = unit(1.25, "cm")))
  dev.off()
