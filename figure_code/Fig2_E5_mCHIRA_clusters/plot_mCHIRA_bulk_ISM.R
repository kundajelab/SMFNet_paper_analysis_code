#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig2_E5_mCHIRA_clusters/plot_mCHIRA_bulk_ISM.R

library(hdf5r)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseqlogo)
library(egg)
library(gridExtra)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(),
                  plot.title = element_text(size = 6), axis.text = element_text(size = 6),
                  title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), plot.margin = grid::unit(c(0, 0, 0, 0), "cm")))

H5_FILE <- "input_files/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_bulk_ISM.h5"

OUT_DIR <- "output/Fig2_E5_mCHIRA_clusters/"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

h5 <- H5File$new(
  H5_FILE,
  mode = "r"
)

regions <- names(h5)

print(regions)

dna_map <- c(
  "0" = "A",
  "1" = "C",
  "2" = "G",
  "3" = "T"
)

influence_list <- list()

for (region in regions) {

  influence <- t(h5[[region]][["influence"]][,])
  dna <- h5[[region]][["variable_dna"]][]

  if (nrow(influence) != length(dna)) {
    influence <- t(influence)
  }

  bulk_influence <- rowSums(
    influence,
    na.rm = TRUE
  )

  seq_letters <- dna_map[
    as.character(dna)
  ]

  influence_list[[region]] <- data.frame(
    region = region,
    pos = seq_along(seq_letters),
    base = seq_letters,
    influence = bulk_influence
  )
}

influence_df <- bind_rows(
  influence_list
)

for (region_use in regions) {
  print(paste0("Plotting ", region_use))

  logo_df <- influence_df %>%
    filter(region == region_use)

  logo_mat <- matrix(
    0,
    nrow = 4,
    ncol = nrow(logo_df),
    dimnames = list(
      c("A", "C", "G", "T"),
      logo_df$pos
    )
  )

  for (i in seq_len(nrow(logo_df))) {

    logo_mat[
      logo_df$base[i],
      i
    ] <- logo_df$influence[i]

  }

  p <- ggseqlogo(
    logo_mat,
    method = "custom",
    seq_type = "dna"
  ) +
    labs(
      title = region_use,
      x = "Position region",
      y = "ISM score"
    ) +
    scale_x_continuous(
      breaks = c(1, 125, 250)
    ) +
    theme_classic() +
    theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          axis.line = element_blank(),
          plot.title = element_text(size = 6), axis.text = element_text(size = 6),
          title = element_text(size = 6),
          axis.title = element_text(size = 6), strip.text = element_text(size = 6),
          strip.background = element_blank(), plot.margin = grid::unit(c(0, 0, 0, 0), "cm"))

    pdf(file.path(
      OUT_DIR,
      paste0("", region_use, "_bulk_influence_logo.pdf")), width = 20/2.54, height = 20/2.54)
  fix <- set_panel_size(p, height = grid::unit(2, "cm"), width = grid::unit(4, "cm"))
  grid.arrange(fix)
  dev.off()
}

write.table(influence_df, file.path(OUT_DIR, "bulk_influence_summed.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
h5$close_all()
