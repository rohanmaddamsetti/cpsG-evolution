#!/usr/bin/env Rscript
# Visualize the manB-cpsG protein alignment using ggmsa.
# Input:  ../results/manB-cpsG_sequences_protein_aligned.fasta (produced by extract_and_align_cpsG_seqs.py)
# Output: ../results/manB-cpsG_protein_alignment.pdf

library(ggmsa)
library(ggplot2)
library(Biostrings)

ALN_FILE    <- "../results/manB-cpsG_sequences_protein_aligned.fasta"
OUTPUT_FILE <- "../results/manB-cpsG_protein_alignment.pdf"

aln     <- readAAMultipleAlignment(ALN_FILE)
aln_len <- ncol(aln)

WRAP_WIDTH <- 60  # positions per panel row

# Rename sequences and reverse order so cpsG-manB chimera appears at top.
aln_ss        <- as(aln, "AAStringSet")
names(aln_ss) <- c("manB", "cpsG", "cpsG-manB chimera")
aln_ss        <- rev(aln_ss)
tmp_file      <- tempfile(fileext = ".fasta")
writeXStringSet(aln_ss, tmp_file)

n_panels   <- ceiling(aln_len / WRAP_WIDTH)
fig_height <- n_panels * 1.3  # ~1.3 inches per panel (3 sequences + x-axis labels)

p <- ggmsa(tmp_file,
           start      = 1,
           end        = aln_len,
           color      = "Chemistry_AA",
           font       = "DroidSansMono",
           char_width = 0.9,
           seq_name   = TRUE,
           border     = "white") +
  facet_msa(field = WRAP_WIDTH) +
  theme(
    axis.text.x      = element_text(size = 8, family = "mono"),
    axis.text.y      = element_text(size = 9, family = "mono", face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.border     = element_rect(fill = NA, color = "grey80", linewidth = 0.4),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    plot.margin      = margin(4, 8, 4, 4),
  )

ggsave(OUTPUT_FILE, p, width = 12, height = fig_height, limitsize = FALSE)
cat("Saved", OUTPUT_FILE, "\n")
