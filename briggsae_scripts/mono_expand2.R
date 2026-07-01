library(rstudioapi)
library(readxl)
library(data.table)
library(openxlsx)

setwd(dirname(getActiveDocumentContext()$path))
cat("Working directory:", getwd(), "\n")

# =============================================================================
# Classify mutations as mono/nonmono using direct perfect mononucleotide repeats
# with +/- 2 bp flanking expansion around each repeat tract.
#
# BED input is standard BED-like:
#   CHROM, start0, end1, repeat_unit
# where start0 is 0-based and end1 is the 1-based half-open BED end.
#
# Expanded interval:
#   start0_exp = max(0, start0 - 2)
#   end1_exp   = end1 + 2
# For 1-based mutation POS matching, this becomes:
#   start1_exp = start0_exp + 1
#   end1_exp   = end1 + 2
#
# The final mono_repeat column below means:
#   TRUE  = POS falls inside repeat tract OR within 2 bp of either end
#   FALSE = POS is outside the expanded repeat interval
# =============================================================================

EXPAND_BP <- 2L

# Load mono repeat BED
bed <- fread("briggsae_direct_perfect_mono_len5.bed",
             col.names = c("CHROM", "start0", "end1", "repeat_unit"),
             header = FALSE)

bed[, repeat_unit := toupper(trimws(repeat_unit))]

# Keep original/core coordinates for optional audit columns
bed[, start0_core := start0]
bed[, end1_core   := end1]

# Expanded coordinates: include 2 bp before start and 2 bp after end
bed[, start0_exp := pmax(0L, start0 - EXPAND_BP)]
bed[, end1_exp   := end1 + EXPAND_BP]
bed[, start1_exp := start0_exp + 1L]  # convert expanded 0-based start to 1-based inclusive

# Core 1-based interval, useful for checking original repeat-only status
bed[, start1_core := start0_core + 1L]

# Separate keyed tables for expanded and core overlap checks
bed_exp <- bed[, .(CHROM, start1_exp, end1_exp, repeat_unit)]
setkey(bed_exp, CHROM, start1_exp, end1_exp)

bed_core <- bed[, .(CHROM, start1_core, end1_core, repeat_unit)]
setkey(bed_core, CHROM, start1_core, end1_core)

classify_mono <- function(mutations) {
  dt <- as.data.table(mutations)
  dt[, row_id := .I]
  dt[, pos_end := POS]  # point interval: POS to POS, 1-based

  # Required for foverlaps
  setkey(dt, CHROM, POS, pos_end)

  # Expanded repeat classification: repeat tract +/- 2 bp
  hits_exp <- foverlaps(
    dt,
    bed_exp,
    by.x = c("CHROM", "POS", "pos_end"),
    by.y = c("CHROM", "start1_exp", "end1_exp"),
    type = "within",
    nomatch = NA
  )

  # Core/original repeat classification: repeat tract only, no flanks
  hits_core <- foverlaps(
    dt,
    bed_core,
    by.x = c("CHROM", "POS", "pos_end"),
    by.y = c("CHROM", "start1_core", "end1_core"),
    type = "within",
    nomatch = NA
  )

  # Use row_id, not row order, because expanded intervals can overlap and
  # foverlaps can return multiple rows for a single mutation.
  exp_ids  <- unique(hits_exp$row_id[!is.na(hits_exp$start1_exp)])
  core_ids <- unique(hits_core$row_id[!is.na(hits_core$start1_core)])

  dt[, mono_repeat_core := row_id %in% core_ids]       # original tract only
  dt[, mono_repeat      := row_id %in% exp_ids]        # tract +/- 2 bp; primary classification
  dt[, mono_flank_only  := mono_repeat & !mono_repeat_core]

  # Optional: annotate repeat unit for expanded hits. If multiple expanded
  # intervals overlap the same mutation, collapse repeat units with ';'.
  unit_map <- hits_exp[!is.na(start1_exp),
                       .(repeat_unit_expanded = paste(sort(unique(repeat_unit)), collapse = ";")),
                       by = row_id]
  dt <- merge(dt, unit_map, by = "row_id", all.x = TRUE, sort = FALSE)

  # Restore original row order and remove helper columns
  setorder(dt, row_id)
  dt[, c("row_id", "pos_end") := NULL]

  dt
}

# Process both files, all sheets
sheets <- c("3x", "5x", "7x", "10x")
files  <- c("SR_final_mutation_list.xlsx", "SR_LR_final_mutation_list.xlsx")

for (f in files) {
  out_name <- sub("\\.xlsx$", "_mono_plus2bp_classified.xlsx", f)
  wb_out   <- openxlsx::createWorkbook()

  for (s in sheets) {
    cat("Classifying", f, "sheet", s, "...\n")
    muts <- read_excel(f, sheet = s)
    classified <- classify_mono(muts)

    cat("  rows:", nrow(classified),
        "| core mono:", sum(classified$mono_repeat_core, na.rm = TRUE),
        "| expanded mono:", sum(classified$mono_repeat, na.rm = TRUE),
        "| flank-only:", sum(classified$mono_flank_only, na.rm = TRUE), "\n")

    openxlsx::addWorksheet(wb_out, s)
    openxlsx::writeData(wb_out, s, classified)
  }

  openxlsx::saveWorkbook(wb_out, out_name, overwrite = TRUE)
  cat("Saved:", out_name, "\n")
}
