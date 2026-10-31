#Run from the repository root: Rscript figure_code/Fig2_E4_PRISM_cluster_benchmarking/plot_CTCF_PRISM.R
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseqlogo)
library(egg)
library(gridExtra)
library(Biostrings)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6),
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

PLOT_DIR <- "output/Fig2_E4_PRISM_cluster_benchmarking/CTCF_PRISM/"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

region <- "chr1_191287950_191289998"
offset <- 874
fasta_file <- "input_files/Fig2_E4_PRISM_cluster_benchmarking/CTCF_example_sequence.fa"

#Raw per-read mutation scores: no scaling or across-read subtraction.
long_df <- read.delim("input_files/Fig2_E4_PRISM_cluster_benchmarking/CTCF_PRISM_scores.tsv") %>%
  mutate(position = as.integer(position_id)) %>%
  arrange(read_id, position)

#Read sequence into R for logoplot
fa <- readDNAStringSet(fasta_file)

#Extract sequence for the region
seq <- fa[[region]]

#Take 300bp starting at offset
seq_window <- subseq(
  seq,
  start = offset + 1,  #Biostrings uses 1-based coordinates
  width = 300
)

#Convert to bases
bases <- strsplit(as.character(seq_window), "")[[1]]

read_list <- split(long_df, long_df$read_id)

#Safety check
pwm_list <- lapply(read_list, function(df) {

  scores <- df$ism_sum

  stopifnot(length(scores) == length(bases))

  pwm <- matrix(0, nrow = 4, ncol = length(bases))
  rownames(pwm) <- c("A", "C", "G", "T")

  for (i in seq_along(bases)) {

    b <- bases[i]
    s <- scores[i]

    if (is.na(s)) next

    pwm[b, i] <- s
  }

  pwm
})

#Plot
all_scores <- long_df$ism_sum
min_val <- min(all_scores, na.rm = TRUE)
max_val <- max(all_scores, na.rm = TRUE)

expand_to_dense <- function(df, value_col = "smf_value", max_pos = 300) {

  df <- df %>%
    arrange(position)

  pos <- df$position
  val <- df[[value_col]]

  mid <- (head(pos, -1) + tail(pos, -1)) / 2

  idx <- findInterval(1:max_pos, c(-Inf, mid, Inf))

  expanded <- data.frame(
    position = 1:max_pos,
    value = val[idx]
  )

  #Run-length encoding
  r <- rle(expanded$value)

  expanded$stretch_id <- rep(seq_along(r$lengths), r$lengths)
  expanded$stretch_length <- rep(r$lengths, r$lengths)

  #Stretch boundaries
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)

  expanded$stretch_start <- rep(starts, r$lengths)
  expanded$stretch_end   <- rep(ends, r$lengths)

  expanded
}

for (cl in names(pwm_list)) {

  #The existing methylation and inferred-state panels use this same read.
  p_logo <- ggseqlogo(pwm_list[[cl]], method = "custom") +
    ggtitle(cl) +
    theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 6),
          axis.ticks.y = element_line(), title = element_blank()) +
    scale_y_continuous(name = "Raw PRISM score", limits = range(c(0, min_val, max_val)))

  expanded_df <- long_df %>%
    filter(read_id == cl) %>%
    filter(smf_value != -1) %>%
    expand_to_dense(value_col = "smf_value", max_pos = 300) %>%
    mutate(edge_val = ifelse(stretch_start == 1 | stretch_end == 300, TRUE, FALSE)) %>%
    mutate(state_id = case_when(
      value == 0 & (stretch_length >= 80 | edge_val == TRUE) ~ 0,
      value == 0 & (stretch_length <= 80 & edge_val == FALSE) ~ 1,
      value == 1 & (stretch_length <= 70 & edge_val == FALSE) ~ 2,
      value == 1 & (stretch_length >= 70 | edge_val == TRUE) ~ 3
    ))

  seg <- expanded_df %>%
    distinct(stretch_id, stretch_start, stretch_end, stretch_length, state_id) %>%
    arrange(stretch_start)

  threshold <- 80

  seg$new_state <- seg$state_id

  seg <- seg %>% arrange(stretch_start)

  allowed <- c(0, 2)

  seg$new_state <- seg$state_id

  for (i in seq_len(nrow(seg))) {

    if (seg$state_id[i] != 1) next

    total_len <- seg$stretch_length[i]

    #Expand left
    j <- i - 1
    while (j >= 1 && seg$state_id[j] %in% allowed) {
      total_len <- total_len + seg$stretch_length[j]
      j <- j - 1
    }

    #Expand right
    j <- i + 1
    while (j <= nrow(seg) && seg$state_id[j] %in% allowed) {
      total_len <- total_len + seg$stretch_length[j]
      j <- j + 1
    }

    #Apply rule
    if (total_len >= 80) {
      seg$new_state[i] <- 0
    }
  }

  p_states <- expanded_df %>%
    left_join(seg) %>%
    ggplot(aes(x = position, fill = as.factor(new_state), y = 1)) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "#e0dede", "1" = "#adadad", "2" = "#159246", "3" = "#000000")) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
          legend.position = "none")

  p_smf <- expanded_df %>%
    left_join(seg) %>%
    ggplot(aes(x = position, y = value)) +
    geom_line() +
    scale_y_continuous(name = "SMF signal", breaks = c(0, 1), labels = c("free", "occupied"))

  pdf(paste0(PLOT_DIR, "CTCF_PRISM_read_", cl, ".pdf"),
      width = 20/2.54, height = 20/2.54)

  p_logo_fixed <- set_panel_size(p_logo, height = unit(1.5, "cm"),width  = unit(4, "cm"))

  p_smf_fixed <- set_panel_size(p_smf, height = unit(0.5, "cm"), width  = unit(4, "cm"))

  p_states_fixed <- set_panel_size(p_states, height = unit(0.25, "cm"), width  = unit(4, "cm"))
  grid.arrange(p_logo_fixed, p_smf_fixed, p_states_fixed, ncol = 1)
  dev.off()

}
