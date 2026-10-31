#Run from SMFNet_paper_analysis_code github repo: Rscript ./figure_code/Fig3_E6_singleMotif_insilico/plot_chip_gc.R

INPUT_DIR <- "input_files/Fig3_E6_singleMotif_insilico"
INTERMEDIATE_DIR <- "output"
RESULTS_DIR <- "output"
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(readr)
library(egg)
library(gridExtra)
library(GenomicRanges)

theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                  axis.line = element_blank(), axis.text.x = element_text(size = 6),
                  title = element_text(size = 6), axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

chip_list <- c(file.path(INPUT_DIR, "ESC_KLF4_PMID28111071.narrowPeak"),
               file.path(INPUT_DIR, "ESC_MYC_PMID25405324.narrowPeak"),
               file.path(INPUT_DIR, "ESC_OCT4_PMID28111071.narrowPeak"),
               file.path(INPUT_DIR, "ESC_SOX2_PMID28111071.narrowPeak"),
               file.path(INPUT_DIR, "ESC_ESRRB_PMID28111071.narrowPeak"),
               file.path(INPUT_DIR, "ESC_REST_PMID22297846.narrowPeak"),
               file.path(INPUT_DIR, "ESC_CTCF_PMID22170606.narrowPeak"),
               file.path(INPUT_DIR, "ESC_ZIC3_PMID20872845.narrowPeak"),
               file.path(INPUT_DIR, "ESC_NFY_PMID25132174.narrowPeak"),
               file.path(INPUT_DIR, "ESC_NRF1_PMID26675734.narrowPeak"),
               file.path(INPUT_DIR, "ESC_BANP_PMID34234345.narrowPeak")
)

PLOT_DIR <- paste0(file.path(RESULTS_DIR, "Fig3_E6_singleMotif_insilico"), "/")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
narrowpeak_cols <- c(
  "chrom", "start", "end", "name", "score",
  "strand", "signalValue", "pValue", "qValue", "peak"
)

read_narrowpeak <- function(path) {

  np <- read_tsv(
    path,
    col_names = narrowpeak_cols,
    col_types = cols(
      chrom = col_character(),
      start = col_double(),
      end = col_double(),
      name = col_character(),
      score = col_double(),
      strand = col_character(),
      signalValue = col_double(),
      pValue = col_double(),
      qValue = col_double(),
      peak = col_double()
    ),
    show_col_types = FALSE
  )

  GRanges(
    seqnames = np$chrom,
    ranges = IRanges(
      start = as.integer(np$start) + 1L,
      end = as.integer(np$end)
    ),
    strand = "*",
    name = np$name,
    score = np$score,
    signalValue = np$signalValue,
    pValue = np$pValue,
    qValue = np$qValue,
    peak = np$peak
  )
}


calc_gc <- function(seq){
  seq <- toupper(as.character(seq))
  chars <- strsplit(seq, "")[[1]]
  mean(chars %in% c("G","C"))
}

gc_list <- list()

for(path in chip_list){

  message(basename(path))

  peaks <- read_narrowpeak(path)

  valid_chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                "chr12", "chr13", "chr14", "chr15", "chr16",  "chr17", "chr18", "chr19", "chrX", "chrY")
  peaks <- peaks[as.character(seqnames(peaks)) %in% valid_chr]

  summit <- start(peaks) + mcols(peaks)$peak

  summit_gr <- GRanges(
    seqnames = seqnames(peaks),
    ranges = IRanges(start = summit, width = 1)
  )

  summit_gr <- resize(summit_gr, width = 1001, fix = "center")

  summit_gr <- trim(summit_gr)

  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, summit_gr)

  gc <- letterFrequency(seqs, c("G","C"), as.prob = TRUE)
  gc <- rowSums(gc)

  gc_list[[basename(path)]] <- tibble::tibble(
    factor = sub("_PMID[0-9]+[.]narrowPeak$", "", basename(path)),
    gc = gc,
    chr = summit_gr$seqnames,
    range = summit_gr$ranges,
    signal_val = mcols(peaks)$signalValue
  )
}

gc_df <- bind_rows(gc_list) %>%
  mutate(factor = if_else(factor == "ESC_NFY", "ESC_NF-YA", factor)) %>%
  group_by(factor) %>%
  mutate(chip_perc = ntile(signal_val, 5),
         gc_perc = ntile(gc, 5)) %>%
  ungroup()

gc_cols <- colorRampPalette( c("#ece2f0", "#a7bddb", "#68a9cf", "#056551"))(11)
gc_order <- gc_df %>%
  group_by(factor) %>%
  summarize(median_gc = median(gc, na.rm = TRUE)) %>%
  arrange(median_gc) %>%
  pull(factor)

gc_perc <- gc_df %>%
  ggplot(aes(x = gc, y = factor(factor, levels = gc_order), fill = factor(factor, levels = gc_order))) +
  geom_violin() +
  stat_summary(geom = "crossbar", fun = "mean", lwd = 0.25, width = 0.5) +
  geom_vline(aes(xintercept = 0.5), linetype = "dashed") +
  coord_cartesian(xlim = c(0.2, 0.8)) +
  scale_fill_manual(values = gc_cols) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.ticks.y = element_blank())

pdf(paste0(PLOT_DIR, "chip_peaks_gc_percs.pdf"),
    width = 20/2.54, height = 20/2.54)
fix <- set_panel_size(gc_perc, height = grid::unit(4, "cm"), width = grid::unit(3, "cm"))
grid.arrange(fix)
dev.off()
