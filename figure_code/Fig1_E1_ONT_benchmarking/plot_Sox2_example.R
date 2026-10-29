library(ggplot2)
library(egg)
library(gridExtra)
library(hdf5r)

#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_sox2_example.R
input_dir <- file.path("input_files", "Fig1_E1_ONT_benchmarking")
output_dir <- file.path("output", "Fig1_E1_ONT_benchmarking")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
avg_dir <- file.path(output_dir, "avg_plots")
sm_dir <- file.path(output_dir, "sm_plots")
dir.create(avg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sm_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- function(sample, kind) {
  file.path(input_dir, paste0("Example_SOX2_locus_", sample, "_", kind, ".h5"))
}

theme_set(theme_classic() + theme(
  legend.text = element_text(size = 6),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line = element_blank(),
  axis.text = element_text(size = 6),
  plot.title = element_text(size = 6),
  title = element_text(size = 6),
  axis.title = element_text(size = 6),
  strip.text = element_text(size = 6),
  strip.background = element_blank(),
  plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
))

save_panel <- function(plot, path, height, width) {
  pdf(path, width = 20 / 2.54, height = 20 / 2.54,
      family = "sans", useDingbats = FALSE)
  on.exit(dev.off())
  panel <- set_panel_size(plot, height = grid::unit(height, "cm"),
                          width = grid::unit(width, "cm"))
  grid.arrange(panel)
}

#Read smfnet.samples matrices; hdf5r reverses the stored matrix dimensions.
read_samples <- function(path, role) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  stopifnot(h5attr(h5, "schema") == "smfnet.samples")
  ids <- as.character(h5[["region_ids"]][])
  result <- setNames(vector("list", length(ids)), ids)
  for (i in seq_along(ids)) {
    group <- h5[[paste0("regions/", sprintf("%08d", i - 1L))]]
    stopifnot(h5attr(group, "region_id") == ids[i])
    coordinates <- strsplit(ids[i], "[:-]")[[1]]
    values <- t(group[[paste0(role, "/smf")]]$read(drop = FALSE))
    positions <- as.numeric(group[[paste0(role, "/positions")]][])
    stopifnot(ncol(values) == length(positions), all(values %in% c(-1, 0, 1)))
    result[[i]] <- list(
      region = ids[i], chrom = coordinates[1], start = as.numeric(coordinates[2]),
      end = as.numeric(coordinates[3]),
      window_start = as.numeric(h5attr(group, "window_start")),
      window_length = length(group[["dna"]][]),
      positions = positions, reads = values
    )
  }
  result
}

#The supplied state files retain the smf_pos/smf_states layout; 0/1 accessible, 2 tf, 3 nucleosome.
read_states <- function(path, samples) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  result <- setNames(vector("list", length(samples)), names(samples))
  for (region in names(samples)) {
    positions <- as.numeric(h5[[paste0("smf_pos/", region)]][])
    values <- t(h5[[paste0("smf_states/", region)]]$read(drop = FALSE))
    stopifnot(identical(positions, samples[[region]]$positions),
              identical(dim(values), dim(samples[[region]]$reads)),
              all(values %in% c(-1, 0, 1, 2, 3)))
    values[values == 1L] <- 0L
    result[[region]] <- values
  }
  result
}

predicted <- read_samples(input_path("SMFNet", "reads"), "predicted")
observed <- read_samples(input_path("Experiment", "reads"), "observed")
stopifnot(setequal(names(predicted), names(observed)))
state_matrices <- list(
  pred = read_states(input_path("SMFNet", "states"), predicted),
  act = read_states(input_path("Experiment", "states"), observed)
)

#Preserve the original site coordinates and inclusive central-window crop.
crop_sample <- function(sample, states, prediction) {
  start <- prediction$start + prediction$window_start
  end <- start + prediction$window_length
  positions <- sample$start + sample$window_start + sample$positions
  keep <- positions >= start & positions <= end
  stopifnot(any(keep))
  list(region = sample$region, positions = positions[keep],
       reads = sample$reads[, keep, drop = FALSE],
       states = states[, keep, drop = FALSE], start = start, end = end)
}

samples <- list(pred = predicted, act = observed)
cropped <- lapply(names(samples), function(source) {
  setNames(lapply(names(predicted), function(region) {
    crop_sample(samples[[source]][[region]], state_matrices[[source]][[region]],
                predicted[[region]])
  }), names(predicted))
})
names(cropped) <- names(samples)

#Compute means from all sampled reads, excluding missing calls separately at each site.
average_profile <- function(sample, source) {
  calls <- sample$reads
  calls[calls == -1L] <- NA_real_
  means <- colMeans(calls, na.rm = TRUE)
  means[is.nan(means)] <- NA_real_
  data.frame(region = sample$region, genomic_pos = sample$positions,
             avg_signal = means, type = source)
}
avg_long <- do.call(rbind, lapply(names(cropped), function(source) {
  do.call(rbind, lapply(cropped[[source]], average_profile, source = source))
}))

average_plot <- function(data, points = FALSE) {
  plot <- ggplot(data, aes(x = genomic_pos, y = avg_signal, color = type,
                          group = interaction(type, region))) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c(act = "#DEDEDE", pred = "#636bb2")) +
    scale_y_continuous(limits = c(0, 1))
  if (points) plot <- plot + geom_point(size = 0.5)
  plot
}

#Only supplied regions are plotted; lines do not bridge unprovided intervals.
sox2_gene <- annotate("rect", xmin = 34650405, xmax = 34652461,
                      ymin = 0, ymax = 0.05, fill = "#DEDEDE", color = "black")
save_panel(average_plot(avg_long) + sox2_gene,
           file.path(avg_dir, "sox2_example_long_avg.pdf"), 2.5, 10)
save_panel(average_plot(avg_long) + sox2_gene +
             scale_x_continuous(limits = c(34640000, 34655000)),
           file.path(avg_dir, "sox2_example_long_avg_15kb_window.pdf"), 2.5, 10)

plot_window <- function(start, end, label) {
  data <- avg_long[avg_long$genomic_pos >= start & avg_long$genomic_pos <= end, ]
  if (!nrow(data)) return(invisible(NULL))
  plot <- average_plot(data, points = TRUE) +
    ggtitle(paste0("Sox2 example | window ", label))
  path <- file.path(avg_dir, paste0("sox2_example_window_", label, "_",
                                  start, "_", end, ".pdf"))
  save_panel(plot, path, 2, 2.5)
}
window_starts <- seq(min(vapply(cropped$pred, `[[`, numeric(1), "start")),
                     max(vapply(cropped$pred, `[[`, numeric(1), "end")) - 500,
                     by = 250)
for (i in seq_along(window_starts)) {
  plot_window(window_starts[i], window_starts[i] + 500, i)
}

examples <- data.frame(
  region = c("chr3:34641671-34643719", "chr3:34646871-34648919",
             "chr3:34648871-34650919"),
  label = c("re", "ctcf", "tss"),
  anchor = c(34642735, 34647892, 34649886)
)
stopifnot(all(examples$region %in% names(predicted)))

#Keep binary-distance/Ward clustering, then group reads by their state nearest the anchor.
plot_molecules <- function(sample, anchor, source, label) {
  calls <- sample$reads
  anchor_index <- which.min(abs(sample$positions - anchor))
  category <- c("undetermined", "accessible", "accessible", "bound", "nuc")[
    sample$states[, anchor_index] + 2L
  ]
  clustered <- if (nrow(calls) > 1L) {
    hclust(dist(calls, method = "binary"), method = "ward.D2")$order
  } else seq_len(nrow(calls))
  category_order <- c("bound", "accessible", "nuc", "undetermined")
  ordered <- clustered[order(match(category[clustered], category_order),
                              seq_along(clustered))]
  read_order <- integer(nrow(calls))
  read_order[ordered] <- seq_along(ordered)
  data <- data.frame(
    smf_pos = rep(sample$positions, times = nrow(calls)),
    reads = as.vector(t(calls)),
    read_id = rep(seq_len(nrow(calls)), each = ncol(calls)),
    order = rep(read_order, each = ncol(calls)),
    cat = rep(category, each = ncol(calls))
  )
  data <- data[data$cat != "undetermined" & data$reads != -1L, ]
  strip <- unique(data[c("read_id", "order", "cat")])
  start <- min(sample$positions)
  end <- max(sample$positions)
  plot <- ggplot(data, aes(x = smf_pos, y = -order, fill = factor(reads))) +
    geom_tile(width = 5) +
    geom_tile(data = strip, aes(x = end - 50, y = -order, fill = cat),
              width = 30, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(start, end)) +
    scale_fill_manual(values = c("0" = "grey", "1" = "black",
                                accessible = "grey", bound = "#0d9146", nuc = "black")) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = "Genomic position", y = "Single molecules", fill = "State")
  save_panel(plot, file.path(sm_dir, paste0(source, "_sox2_", label,
                                          "_example_sm_reads.pdf")), 1.5, 2.5)
  invisible(list(order = ordered, category = category))
}

for (i in seq_len(nrow(examples))) {
  example <- examples[i, ]
  prediction <- cropped$pred[[example$region]]
  plot_window(prediction$start, prediction$end, paste0("example_", example$label))
  for (source in names(cropped)) {
    plot_molecules(cropped[[source]][[example$region]], example$anchor,
                   source, example$label)
  }
}
