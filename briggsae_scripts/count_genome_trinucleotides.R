# =============================================================================
# Step 3: Count all trinucleotides in the C. briggsae reference genome
#         (replaces briggsae_count_trinucleotides.py)
#
# Run on HiPerGator (needs Biostrings -- available in Bioconductor)
#   install if needed: BiocManager::install("Biostrings")
#
# Input:  20250626_c_briggsae_Feb2020.genome.fa
# Output: briggsae_genome_trinucleotide_counts.csv
#
# This is the normalization denominator for per-context mutation rates.
# Takes ~2-5 minutes for the briggsae genome.
# =============================================================================
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(rstudioapi)
library(Biostrings)
library(dplyr)

setwd(dirname(getActiveDocumentContext()$path))

FASTA <- "20250626_c_briggsae_Feb2020.genome.fa"
OUT   <- "briggsae_genome_trinucleotide_counts.csv"

# ── Load genome ───────────────────────────────────────────────────────────────
cat("Loading reference genome...\n")
genome <- readDNAStringSet(FASTA)
cat(sprintf("  %d chromosomes / scaffolds\n", length(genome)))
cat(sprintf("  Total length: %s bp\n",
            format(sum(width(genome)), big.mark = ",")))

# ── Count trinucleotides per chromosome, sum across genome ───────────────────
# trinucleotideFrequency() from Biostrings counts all 64 3-mers per sequence
# step=1 means overlapping counts (standard for trinucleotide context)

cat("Counting trinucleotides (overlapping, step=1)...\n")
tri_mat <- trinucleotideFrequency(genome, step = 1)   # rows = chroms, cols = 64 contexts
tri_totals <- colSums(tri_mat)                         # sum across all chromosomes

cat(sprintf("  Total trinucleotides counted: %s\n",
            format(sum(tri_totals), big.mark = ",")))

# ── Write output ──────────────────────────────────────────────────────────────
tri_df <- data.frame(
  trinucleotide = names(tri_totals),
  count         = as.integer(tri_totals),
  stringsAsFactors = FALSE
)

write.csv(tri_df, OUT, row.names = FALSE)
cat(sprintf("Saved: %s  (%d trinucleotide contexts)\n", OUT, nrow(tri_df)))
cat("\nReady to run local_context.R\n")
