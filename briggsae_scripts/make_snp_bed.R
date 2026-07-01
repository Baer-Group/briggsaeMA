# =============================================================================
# Step 1: Create mutations_snp_3bp.bed from mutation xlsx
#         (replaces make_snp_bed.py)
#
# Run this on HiPerGator before bedtools getfasta and local_context.R
#
# Inputs:
#   SR_final_mutation_list_mono_classified.xlsx
#
# Outputs:
#   mutations_snp_3bp.bed   -- ready for bedtools getfasta
#
# BED format (0-based, 3-bp window):
#   CHROM   POS-2   POS+1   CHROM_POS_REF_ALT_Sample
#
# After running this script, run on HiPerGator:
#   bedtools getfasta \
#       -fi 20250626_c_briggsae_Feb2020.genome.fa \
#       -bed mutations_snp_3bp.bed \
#       -fo mutations_snp_3bp.fasta \
#       -name
# =============================================================================

library(rstudioapi)
library(readxl)
library(dplyr)

setwd(dirname(getActiveDocumentContext()$path))

SHEET <- "3x"

# ── Load mutation list ────────────────────────────────────────────────────────
cat("Loading mutation list...\n")
muts <- read_excel("SR_final_mutation_list_mono_classified.xlsx", sheet = SHEET)
cat(sprintf("  %d rows loaded\n", nrow(muts)))

# ── Filter to SNVs only ───────────────────────────────────────────────────────
snvs <- muts %>%
  filter(TYPE == "snp") %>%
  filter(!is.na(CHROM), !is.na(POS), !is.na(REF), !is.na(ALT))

cat(sprintf("  %d SNV mutations after filtering\n", nrow(snvs)))

# ── Build BED ─────────────────────────────────────────────────────────────────
# BED is 0-based half-open: [start, end)
# VCF POS is 1-based → 3-bp window:
#   bed_start = POS - 2   (0-based, one base BEFORE the mutation)
#   bed_end   = POS + 1   (exclusive, one base AFTER the mutation)
# Middle base of extracted 3-mer = the mutation site (REF base)

bed <- snvs %>%
  mutate(
    bed_start = as.integer(POS) - 2L,
    bed_end   = as.integer(POS) + 1L,
    # Unique name for joining back in local_context.R
    name      = paste(CHROM, POS, REF, ALT, Sample, sep = "_")
  ) %>%
  filter(bed_start >= 0) %>%          # drop any at chromosome start
  select(CHROM, bed_start, bed_end, name)

cat(sprintf("  %d BED records written (%d dropped: start < 0)\n",
            nrow(bed), nrow(snvs) - nrow(bed)))

# ── Write BED (no header, tab-separated) ─────────────────────────────────────
write.table(bed, "mutations_snp_3bp.bed",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

cat("Saved: mutations_snp_3bp.bed\n")
cat("\nNext step on HiPerGator:\n")
cat("  bedtools getfasta \\\n")
cat("      -fi 20250626_c_briggsae_Feb2020.genome.fa \\\n")
cat("      -bed mutations_snp_3bp.bed \\\n")
cat("      -fo mutations_snp_3bp.fasta \\\n")
cat("      -name\n")
