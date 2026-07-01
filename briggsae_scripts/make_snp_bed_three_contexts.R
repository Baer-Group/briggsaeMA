# =============================================================================
# Create SNP 3-bp context BED files for three orientation checks
#
# Inputs:
#   SR_final_mutation_list_mono_classified.xlsx
#
# Outputs:
#   mutations_snp_3bp_center_xYz.bed       # Rajaei-style: mutation in middle
#   mutations_snp_3bp_3prime_xxY.bed       # Saxena Fig 3B-style: mutation at 3' end
#   mutations_snp_3bp_5prime_Yxx.bed       # complementary check: mutation at 5' end
#
# Coordinate logic:
#   VCF POS is 1-based.
#   BED is 0-based, half-open [start, end).
#
#   center xYz: positions POS-1, POS, POS+1
#       BED start = POS - 2
#       BED end   = POS + 1
#       Expected check after getfasta: substr(context, 2, 2) == REF
#
#   3prime xxY: positions POS-2, POS-1, POS
#       BED start = POS - 3
#       BED end   = POS
#       Expected check after getfasta: substr(context, 3, 3) == REF
#
#   5prime Yxx: positions POS, POS+1, POS+2
#       BED start = POS - 1
#       BED end   = POS + 2
#       Expected check after getfasta: substr(context, 1, 1) == REF
# =============================================================================

library(rstudioapi)
library(readxl)
library(dplyr)

setwd(dirname(getActiveDocumentContext()$path))

SHEET <- "3x"
MUTATION_XLSX <- "SR_final_mutation_list_mono_classified.xlsx"

cat("Loading mutation list...\n")
muts <- read_excel(MUTATION_XLSX, sheet = SHEET)
cat(sprintf("  %d rows loaded\n", nrow(muts)))

snvs <- muts %>%
  filter(TYPE == "snp") %>%
  filter(!is.na(CHROM), !is.na(POS), !is.na(REF), !is.na(ALT), !is.na(Sample)) %>%
  mutate(
    POS = as.integer(POS),
    REF = toupper(as.character(REF)),
    ALT = toupper(as.character(ALT)),
    name = paste(CHROM, POS, REF, ALT, Sample, sep = "_")
  ) %>%
  filter(nchar(REF) == 1, nchar(ALT) == 1,
         REF %in% c("A", "C", "G", "T"),
         ALT %in% c("A", "C", "G", "T"))

cat(sprintf("  %d SNV mutations after filtering\n", nrow(snvs)))

write_bed <- function(df, bed_start, bed_end, outfile) {
  bed <- df %>%
    mutate(
      bed_start = bed_start,
      bed_end   = bed_end
    ) %>%
    filter(bed_start >= 0) %>%
    select(CHROM, bed_start, bed_end, name)

  write.table(
    bed, outfile,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )

  cat(sprintf("  Saved %-38s : %d records (%d dropped because start < 0)\n",
              outfile, nrow(bed), nrow(df) - nrow(bed)))
}

cat("Writing BED files...\n")

# 1. Rajaei-style centered motif: 5'-xYz-3', Y is mutable base
write_bed(
  snvs,
  as.integer(snvs$POS) - 2L,
  as.integer(snvs$POS) + 1L,
  "mutations_snp_3bp_center_xYz.bed"
)

# 2. Saxena-style 3' motif: 5'-xxY-3', Y is mutable base at 3' end
write_bed(
  snvs,
  as.integer(snvs$POS) - 3L,
  as.integer(snvs$POS),
  "mutations_snp_3bp_3prime_xxY.bed"
)

# 3. 5' motif: 5'-Yxx-3', Y is mutable base at 5' end
write_bed(
  snvs,
  as.integer(snvs$POS) - 1L,
  as.integer(snvs$POS) + 2L,
  "mutations_snp_3bp_5prime_Yxx.bed"
)

cat("\nRun these bedtools commands on HiPerGator:\n\n")
cat("GENOME=20250626_c_briggsae_Feb2020.genome.fa\n\n")

cat("bedtools getfasta \\\n")
cat("  -fi $GENOME \\\n")
cat("  -bed mutations_snp_3bp_center_xYz.bed \\\n")
cat("  -fo mutations_snp_3bp_center_xYz.fasta \\\n")
cat("  -name\n\n")

cat("bedtools getfasta \\\n")
cat("  -fi $GENOME \\\n")
cat("  -bed mutations_snp_3bp_3prime_xxY.bed \\\n")
cat("  -fo mutations_snp_3bp_3prime_xxY.fasta \\\n")
cat("  -name\n\n")

cat("bedtools getfasta \\\n")
cat("  -fi $GENOME \\\n")
cat("  -bed mutations_snp_3bp_5prime_Yxx.bed \\\n")
cat("  -fo mutations_snp_3bp_5prime_Yxx.fasta \\\n")
cat("  -name\n\n")

cat("Recommended validation after reading FASTA:\n")
cat("  center_xYz : substr(context, 2, 2) == REF\n")
cat("  3prime_xxY : substr(context, 3, 3) == REF\n")
cat("  5prime_Yxx : substr(context, 1, 1) == REF\n")
cat("Also drop any FASTA records where nchar(context) != 3 or context contains N.\n")

