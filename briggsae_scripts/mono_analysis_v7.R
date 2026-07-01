# =============================================================================
# Mononucleotide Repeat Analysis -- C. briggsae MA Lines
# Framework: Rajaei et al. 2021 (Genome Research)
# v4: context-corrected SNV denominators; outputs use *_ctxcorr filenames
# =============================================================================
#
# KEY METHODOLOGICAL POINTS FROM RAJAEI ET AL. 2021:
#
#   1. Compute mutation RATES separately within mono and non-mono contexts
#      (NOT enrichment ratios vs genome proportion -- that was wrong)
#
#   2. Context-specific callable sites:
#        callable_mono    = callable_total x P_MONO
#        callable_nonmono = callable_total x (1 - P_MONO)
#      where P_MONO = REPEAT_BP / GENOME_SIZE = 6001892 / 106196309
#
#   3. Rate formula per line i, class c, context ctx:
#        rate_i_c_ctx = count_i_c_ctx / (callable_ctx_i x frac_c x generations_i)
#        frac_c = GC_FRAC for GC-context SNVs, AT_FRAC for AT-context SNVs,
#                 1.0 for indels (Rajaei Table 1 normalisation)
#
#   4. Primary output = rate ratio:  mu_mono / mu_nonmono  (Rajaei key result:
#        AT->TA is ~7x greater in mono; +/-1bp indels are ~26-39x greater)
#
#   5. Spectrum comparisons use Monte Carlo Fisher's Exact Test (MC-FET)
#        on 6-class count matrices -- same as Rajaei et al.
#
#   6. Per-unit (A/T/C/G run) analysis is an EXTENSION beyond Rajaei
#
# Inputs:
#   master_table.csv                               -- per-line callable, gens
#   SR_final_mutation_list_mono_classified.xlsx    -- mono_repeat col present
#   SR_LR_final_mutation_list_mono_classified.xlsx
#   briggsae_mono_repeats_with_unit.bed            -- 4-col: CHROM start0 end1 unit
#
# Outputs (Rajaei-equivalent):
#   rate_table_rajaei.csv     -- Table 1 equiv: nonmono / mono / total rates
#   rate_ratio_bootstrap.csv  -- mu_mono/mu_nonmono with 95% CI
#   spectrum_stats.csv        -- MC-FET p-values for spectrum comparisons
#   Fig_SNV_rates_plus2_olddenom_monoATGC_unitAssoc.pdf/.png    -- Fig 2A-C equiv: SNV rates by context
#   Fig_Indel_rates_plus2_olddenom_monoATGC_unitAssoc.pdf/.png  -- Fig 2D-F equiv: indel rates by context
#   Fig_SNV_spectrum_plus2_olddenom_monoATGC_unitAssoc.pdf/.png -- Fig 4A-C equiv: proportions by context
#   Fig_rate_ratio_plus2_olddenom_monoATGC_unitAssoc.pdf/.png   -- summary: rate ratio with bootstrap CI
#   Fig_TsTv_context_plus2_olddenom_monoATGC_unitAssoc.pdf/.png
#   Fig_unit_heatmap_plus2_olddenom_monoATGC_unitAssoc.pdf/.png -- EXTENSION: per-unit enrichment
# =============================================================================

library(rstudioapi)
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd(dirname(getActiveDocumentContext()$path))

# =============================================================================
# CONSTANTS
# =============================================================================
GENOME_SIZE <- 106196309
REPEAT_BP   <- 8360213
P_MONO      <- REPEAT_BP / GENOME_SIZE  
P_NONMONO   <- 1 - P_MONO              

GC_FRAC     <- 0.37
AT_FRAC     <- 1 - GC_FRAC              # 0.63
GENERATIONS <- 238
N_BOOT      <- 10000
set.seed(2025)

SHEET       <- "3x"
OUT_SUFFIX  <- "plus2_olddenom_monoATGC_unitAssoc"

# SNV classes -- consistent with spectra.R and Rajaei Table 1
SNP_CLASSES <- c("GC_to_AT","GC_to_CG","GC_to_TA",
                 "AT_to_GC","AT_to_CG","AT_to_TA")
SNP_LABELS  <- c("GC->AT","GC->CG","GC->TA",
                 "AT->GC","AT->CG","AT->TA")
# GC/AT fraction for rate normalisation (Rajaei Table 1 approach)
SNP_FRAC    <- c(GC_FRAC, GC_FRAC, GC_FRAC,
                 AT_FRAC, AT_FRAC, AT_FRAC)
names(SNP_FRAC) <- SNP_CLASSES

TS_CLASSES  <- c("GC_to_AT","AT_to_GC")   # transitions

# Indel bins -- signed format
INDEL_BINS   <- c("del_10plus","del_6_10","del_3_5","del_2","del_1",
                  "ins_1","ins_2","ins_3_5","ins_6_10","ins_10plus")
INDEL_LABELS <- c("<-10","-10 to -6","-5 to -3","-2","-1",
                  "1","2","3 to 5","6 to 10","10<")

ALL_CLASSES  <- c(SNP_CLASSES, INDEL_BINS)

UNIT_LEVELS  <- c("A","T","C","G")
UNIT_COLORS  <- c(A="#4DAF4A",T="#E41A1C",C="#377EB8",G="#FF7F00")

COL_HK   <- "blue4";  COL_PB   <- "coral"
CONTEXTS <- c("nonmono","mono","total")
CTX_COLORS <- c(nonmono="grey50", mono="#D62728", total="grey20")

# =============================================================================
# SHARED HELPERS
# =============================================================================
sem <- function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x)))

theme_rajaei <- function(base_size=18) {
  theme_bw(base_size=base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x  = element_text(angle=45, hjust=1, size=18),
      axis.text.y  = element_text(size=18),
      axis.title   = element_text(size=20),
      legend.position  = "top",
      legend.title     = element_text(size=18),
      legend.text      = element_text(size=18),
      strip.background = element_rect(fill="grey92"),
      strip.text       = element_text(size=18, face="bold"),
      plot.subtitle    = element_text(size=20, face="bold"),
      plot.title       = element_blank()
    )
}

save_fig <- function(p, name, w, h) {
  ggsave(paste0(name,".pdf"), p, width=w, height=h)
  ggsave(paste0(name,".png"), p, width=w, height=h, dpi=300)
  cat(sprintf("  Saved: %s.pdf/.png\n", name))
}

# =============================================================================
# STEP 1: LOAD MASTER TABLE
# Source of callable sites and generations per line
# =============================================================================
cat("Loading master table...\n")

mt <- read.csv("master_table.csv", stringsAsFactors=FALSE)

# Strain assignment (just in case Strain col is absent)
if (!"Strain" %in% names(mt)) {
  mt$sample_num <- as.integer(sub("_.*","", as.character(mt$Line)))
  mt$Strain     <- ifelse(mt$sample_num <= 298, "HK104", "PB800")
}

mt <- mt %>%
  filter(!is.na(SR_callable_3x)) %>%
  select(Line, Strain, callable = SR_callable_3x) %>%
  mutate(Line = as.character(Line))

cat(sprintf("  %d lines: HK104=%d, PB800=%d\n",
            nrow(mt),
            sum(mt$Strain=="HK104"),
            sum(mt$Strain=="PB800")))

# =============================================================================
# STEP 2: LOAD 4-COLUMN BED -- per-unit null for extension
# =============================================================================
cat("Loading repeat BED file...\n")

bed <- fread("briggsae_direct_perfect_mono_len5.bed",
             col.names=c("CHROM","start0","end1","repeat_unit"),
             header=FALSE)
bed[, repeat_unit := toupper(trimws(repeat_unit))]
bed[, start1 := start0 + 1]
setkey(bed, CHROM, start1, end1)

# Expanded intervals used ONLY to assign repeat_unit to mutations that fall
# within the original repeat tract +/- 2 bp. The mono/nonmono status itself
# is read from the plus2-classified Excel files.
bed_exp <- copy(bed)
bed_exp[, start1_exp := pmax(1L, start1 - 2L)]
bed_exp[, end1_exp   := end1 + 2L]
setkey(bed_exp, CHROM, start1_exp, end1_exp)

# Per-unit bp for denominator/composition analysis uses the original core repeat tracts.
unit_bp <- bed[, .(unit_bp = sum(end1-start0)), by=repeat_unit]
unit_bp[, p_null_unit := unit_bp / GENOME_SIZE]
cat("Per-unit null from core repeat tracts:\n"); print(unit_bp)

# Per-unit bp for the heatmap null uses repeat-associated intervals: core tract +/- 2 bp.
# This is only an in-memory object; it does not modify the BED file on disk.
# It makes the heatmap definition match the mono_repeat numerator used in the barplot.
unit_bp_assoc <- bed_exp[, .(unit_bp_assoc = sum(end1_exp - start1_exp + 1L)), by=repeat_unit]
unit_bp_assoc[, p_null_assoc := unit_bp_assoc / sum(unit_bp_assoc)]
cat("Per-unit null for repeat-associated heatmap, tract +/- 2 bp:\n"); print(unit_bp_assoc)

# -----------------------------------------------------------------------------
# Context-corrected nucleotide denominators for SNV rates
# -----------------------------------------------------------------------------
# The mono-repeat compartment is not GC/AT-balanced. Therefore, SNV denominators
# must use GC/AT bp fractions separately within mono and nonmono compartments,
# not P_MONO * genome-wide GC_FRAC/AT_FRAC.
get_unit_bp <- function(u) {
  val <- unit_bp$unit_bp[unit_bp$repeat_unit == u]
  if (length(val) == 0) 0 else as.numeric(val[1])
}

MONO_A_BP  <- get_unit_bp("A")
MONO_T_BP  <- get_unit_bp("T")
MONO_C_BP  <- get_unit_bp("C")
MONO_G_BP  <- get_unit_bp("G")
MONO_AT_BP <- MONO_A_BP + MONO_T_BP
MONO_GC_BP <- MONO_C_BP + MONO_G_BP
MONO_AT_PROP <- MONO_AT_BP / (MONO_AT_BP + MONO_GC_BP)
MONO_GC_PROP <- MONO_GC_BP / (MONO_AT_BP + MONO_GC_BP)

# Old fixed mono opportunity, split by AT/GC composition inside mono repeats.
# This is the compromise you requested:
#   classification/numerator = original tract +/- 2 bp from Excel mono_repeat
#   total mono denominator   = old fixed REPEAT_BP / GENOME_SIZE
#   SNV mono denominator     = P_MONO split by mono AT% and mono GC%
MONO_AT_FRAC <- P_MONO * MONO_AT_PROP
MONO_GC_FRAC <- P_MONO * MONO_GC_PROP
NONMONO_AT_FRAC <- AT_FRAC - MONO_AT_FRAC
NONMONO_GC_FRAC <- GC_FRAC - MONO_GC_FRAC

TOTAL_GC_BP <- GENOME_SIZE * GC_FRAC
TOTAL_AT_BP <- GENOME_SIZE * AT_FRAC

if (NONMONO_GC_FRAC <= 0 || NONMONO_AT_FRAC <= 0) {
  stop("Invalid denominator fractions. Check REPEAT_BP, GC_FRAC/AT_FRAC, or repeat BED.")
}

DENOM_FRAC <- tibble::tibble(
  base_context = c("GC", "AT", "indel"),
  mono_frac    = c(MONO_GC_FRAC, MONO_AT_FRAC, P_MONO),
  nonmono_frac = c(NONMONO_GC_FRAC, NONMONO_AT_FRAC, P_NONMONO),
  total_frac   = c(GC_FRAC, AT_FRAC, 1.0)
)

cat("\nContext-corrected denominator fractions:\n")
print(DENOM_FRAC)

# =============================================================================
# STEP 3: LOAD MUTATION LISTS + ADD repeat_unit via foverlaps
# =============================================================================
cat("\nLoading and annotating mutation lists...\n")

PYRIM_MAP <- c(
  C_to_T="GC_to_AT", C_to_G="GC_to_CG", C_to_A="GC_to_TA",
  T_to_C="AT_to_GC", T_to_G="AT_to_CG", T_to_A="AT_to_TA"
)

classify_snv <- function(ref, alt) {
  comp <- c(A="T",T="A",C="G",G="C")
  r <- toupper(ref); a <- toupper(alt)
  if (r %in% c("G","A")) { r <- comp[r]; a <- comp[a] }
  PYRIM_MAP[paste0(r,"_to_",a)]
}

classify_indel <- function(ref, alt) {
  d <- nchar(alt) - nchar(ref)
  if      (d <= -10) "del_10plus"
  else if (d <=  -6) "del_6_10"
  else if (d <=  -3) "del_3_5"
  else if (d ==  -2) "del_2"
  else if (d ==  -1) "del_1"
  else if (d ==   1) "ins_1"
  else if (d ==   2) "ins_2"
  else if (d <=   5) "ins_3_5"
  else if (d <=  10) "ins_6_10"
  else               "ins_10plus"
}

# Robust TRUE/FALSE conversion. This fixes warnings/NA values when Excel stores
# logical values as TRUE/FALSE text rather than 1/0.
to_logical_clean <- function(x) {
  if (is.logical(x)) return(x)
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(z))
  out[z %in% c("true","t","1","yes","y")]  <- TRUE
  out[z %in% c("false","f","0","no","n")] <- FALSE
  out
}

find_input_file <- function(primary, fallback=NULL) {
  candidates <- c(primary, fallback)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  stop(
    "Cannot find input Excel file. Looked for: ",
    paste(candidates, collapse=" ; "),
    "\nCurrent working directory: ", getwd(),
    "\nExcel files here: ", paste(list.files(pattern="\\.xlsx$"), collapse=" ; ")
  )
}

load_mutations <- function(file) {
  file <- normalizePath(file, mustWork=TRUE)
  df <- as.data.table(read_excel(file, sheet=SHEET))

  # Strain
  df[, sample_num := suppressWarnings(as.integer(sub("_.*","", as.character(Sample))))]
  df[, Strain := ifelse(sample_num <= 298, "HK104", "PB800")]
  df[, Line   := as.character(sample_num)]
  df[, sample_num := NULL]

  # foverlaps for repeat_unit using original tract +/- 2 bp.
  # This makes the heatmap/unit analyses use the same plus2 context as the
  # mono_repeat column. If a mutation overlaps multiple expanded intervals,
  # keep the first unit only for plotting; mono/nonmono status is unaffected.
  df[, row_id := .I]
  df[, pos_end := POS]
  q <- df[, .(row_id, CHROM, POS, pos_end)]
  setkey(q, CHROM, POS, pos_end)

  hits <- foverlaps(q,
                    bed_exp[, .(CHROM, start1_exp, end1_exp, repeat_unit, core_start1=start1, core_end1=end1)],
                    by.x=c("CHROM","POS","pos_end"),
                    by.y=c("CHROM","start1_exp","end1_exp"),
                    type="within", nomatch=NA)

  # Assign exactly one repeat-associated unit to each mutation.
  # If multiple expanded tracts overlap one mutation, prefer a true core hit;
  # otherwise choose the nearest core tract. This prevents double-counting.
  hits[, in_core := !is.na(core_start1) & POS >= core_start1 & POS <= core_end1]
  hits[, dist_to_core := fifelse(
    is.na(core_start1), Inf,
    fifelse(POS < core_start1, core_start1 - POS,
            fifelse(POS > core_end1, POS - core_end1, 0))
  )]

  hit_ids <- unique(hits[!is.na(repeat_unit), row_id])
  unit_map <- hits[!is.na(repeat_unit)][
    order(row_id, -in_core, dist_to_core),
    .SD[1], by=row_id
  ][, .(row_id, repeat_unit, repeat_context_source = ifelse(in_core, "core", "plus2_flank"))]

  df[, repeat_unit := NA_character_]
  df[, repeat_context_source := NA_character_]
  df[unit_map, `:=`(repeat_unit = i.repeat_unit,
                    repeat_context_source = i.repeat_context_source), on="row_id"]

  # This is only a diagnostic check. We still use the Excel mono_repeat column.
  if ("mono_repeat" %in% names(df)) {
    mono_file <- to_logical_clean(df$mono_repeat)
    n_discord <- sum(!is.na(mono_file) & (mono_file != (df$row_id %in% hit_ids)))
    if (n_discord > 0) {
      warning(sprintf("%s: %d rows differ between Excel mono_repeat and direct plus2 overlap check.",
                      basename(file), n_discord))
    }
  }

  df[, c("row_id","pos_end") := NULL]

  # Classify mutation class
  df <- df %>%
    rowwise() %>%
    mutate(mut_class = case_when(
      TYPE=="snp"   ~ classify_snv(REF, ALT),
      TYPE=="indel" ~ classify_indel(REF, ALT),
      TRUE          ~ NA_character_
    )) %>%
    ungroup()

  # Ensure mono_repeat is logical. This uses the plus2 classification from Excel.
  if (!"mono_repeat" %in% names(df)) {
    stop(basename(file), " has no mono_repeat column. Use the plus2-classified Excel file.")
  }
  df$mono_repeat <- to_logical_clean(df$mono_repeat)
  if ("mono_repeat_core" %in% names(df)) df$mono_repeat_core <- to_logical_clean(df$mono_repeat_core)
  if ("mono_flank_only"  %in% names(df)) df$mono_flank_only  <- to_logical_clean(df$mono_flank_only)
  df
}

SR_FILE <- find_input_file("SR_final_mutation_list_mono_plus2bp_classified.xlsx",
                           "SR_final_mutation_list_mono_plus2bp_RECLASSIFIED.xlsx")
SR_LR_FILE <- find_input_file("SR_LR_final_mutation_list_mono_plus2bp_classified.xlsx",
                              "SR_LR_final_mutation_list_mono_plus2bp_RECLASSIFIED.xlsx")

cat("Using SR file:    ", SR_FILE, "\n", sep="")
cat("Using SR+LR file: ", SR_LR_FILE, "\n", sep="")

sr    <- load_mutations(SR_FILE)
sr_lr <- load_mutations(SR_LR_FILE)

cat(sprintf("  SR:    %d mutations | %d in plus2 mono context\n",
            nrow(sr), sum(sr$mono_repeat, na.rm=TRUE)))
cat(sprintf("  SR+LR: %d mutations | %d in plus2 mono context\n",
            nrow(sr_lr), sum(sr_lr$mono_repeat, na.rm=TRUE)))

cat("\nSR 3x diagnostic counts by TYPE and plus2 mono_repeat:\n")
print(sr %>% count(TYPE, mono_repeat))
cat("\nSR+LR 3x diagnostic counts by TYPE and plus2 mono_repeat:\n")
print(sr_lr %>% count(TYPE, mono_repeat))

datasets <- list(SR=sr, SR_LR=sr_lr)

# =============================================================================
# STEP 4: AGGREGATE PER-LINE x CLASS x CONTEXT COUNTS
#
# For each line, for each mutation class, count mutations in:
#   mono    (mono_repeat == TRUE)
#   nonmono (mono_repeat == FALSE)
#   total   (both combined)
#
# Lines with 0 mutations of a class are included via complete grid expansion.
# This is essential -- zero-count lines contribute to the mean rate.
# =============================================================================
cat("\nAggregating per-line x class x context counts...\n")

aggregate_counts <- function(mut_df, mt_df) {

  # Count per line x class x context from mutation list
  raw_counts <- mut_df %>%
    filter(!is.na(mut_class), !is.na(mono_repeat)) %>%
    mutate(context = ifelse(mono_repeat, "mono", "nonmono")) %>%
    group_by(Line, Strain, mut_class, context) %>%
    summarise(count = n(), .groups="drop")

  # Full grid: every line x every class x both contexts
  # Lines not present for a class x context = 0 mutations
  all_lines  <- mt_df %>% select(Line, Strain)
  all_combos <- expand.grid(
    Line      = all_lines$Line,
    mut_class = ALL_CLASSES,
    context   = c("mono","nonmono"),
    stringsAsFactors = FALSE
  ) %>%
    left_join(all_lines, by="Line") %>%
    left_join(raw_counts, by=c("Line","Strain","mut_class","context")) %>%
    replace_na(list(count=0))

  # Add total context (mono + nonmono combined)
  total_counts <- all_combos %>%
    group_by(Line, Strain, mut_class) %>%
    summarise(count=sum(count), .groups="drop") %>%
    mutate(context="total")

  bind_rows(all_combos, total_counts) %>%
    left_join(mt_df %>% select(Line, callable), by="Line")
}

counts_sr    <- aggregate_counts(sr,    mt)
counts_sr_lr <- aggregate_counts(sr_lr, mt)

cat(sprintf("  SR grid:    %d rows\n", nrow(counts_sr)))
cat(sprintf("  SR+LR grid: %d rows\n", nrow(counts_sr_lr)))

# =============================================================================
# STEP 5: COMPUTE PER-LINE RATES (Rajaei Table 1 formula)
#
# Rate formula (same as Rajaei et al. 2021):
#   SNV: rate = count / (callable_ctx x GC_or_AT_frac x generations)
#   Indel: rate = count / (callable_ctx x generations)     [no GC/AT norm]
#
# callable_ctx:
#   mono    context: callable x P_MONO
#   nonmono context: callable x P_NONMONO
#   total   context: callable
#
# P_MONO = 6001892 / 106196309 = 0.05652
# =============================================================================
cat("Computing per-line rates...\n")

compute_rates <- function(counts_df) {
  counts_df %>%
    mutate(
      # Base context for denominator choice
      base_context = case_when(
        mut_class %in% c("GC_to_AT", "GC_to_CG", "GC_to_TA") ~ "GC",
        mut_class %in% c("AT_to_GC", "AT_to_CG", "AT_to_TA") ~ "AT",
        TRUE ~ "indel"
      ),
      # Correct context-specific callable fraction:
      #   GC-class mono SNVs use mono GC-repeat bp / genome,
      #   AT-class mono SNVs use mono AT-repeat bp / genome,
      #   nonmono SNVs use the corresponding remainder,
      #   indels use total mono/nonmono bp fraction.
      denom_frac = case_when(
        context == "mono"    & base_context == "GC"    ~ MONO_GC_FRAC,
        context == "mono"    & base_context == "AT"    ~ MONO_AT_FRAC,
        context == "mono"    & base_context == "indel" ~ P_MONO,
        context == "nonmono" & base_context == "GC"    ~ NONMONO_GC_FRAC,
        context == "nonmono" & base_context == "AT"    ~ NONMONO_AT_FRAC,
        context == "nonmono" & base_context == "indel" ~ P_NONMONO,
        context == "total"   & base_context == "GC"    ~ GC_FRAC,
        context == "total"   & base_context == "AT"    ~ AT_FRAC,
        context == "total"   & base_context == "indel" ~ 1.0,
        TRUE ~ NA_real_
      ),
      callable_ctx = callable * denom_frac,
      rate = count / (callable_ctx * GENERATIONS)
    )
}
rates_sr    <- compute_rates(counts_sr)
rates_sr_lr <- compute_rates(counts_sr_lr)

# =============================================================================
# STEP 6: SUMMARISE RATES (mean +/- SEM per Strain x class x context)
# This is the Rajaei Table 1 equivalent
# =============================================================================
cat("Summarising rates...\n")

summarise_rates <- function(rates_df, ds_label) {
  rates_df %>%
    group_by(Strain, mut_class, context) %>%
    summarise(
      n_lines   = n(),
      mean_rate = mean(rate, na.rm=TRUE),
      se_rate   = sem(rate),
      .groups   = "drop"
    ) %>%
    mutate(dataset = ds_label)
}

rate_summary_sr    <- summarise_rates(rates_sr,    "SR")
rate_summary_sr_lr <- summarise_rates(rates_sr_lr, "SR_LR")
all_rate_summary   <- bind_rows(rate_summary_sr, rate_summary_sr_lr)

# Wide format matching Rajaei Table 1:
# Strain | Class | Nonmono (SEM) | Mono (SEM) | Total (SEM)
rate_table <- all_rate_summary %>%
  filter(dataset=="SR") %>%
  select(Strain, mut_class, context, mean_rate, se_rate) %>%
  pivot_wider(
    names_from  = context,
    values_from = c(mean_rate, se_rate)
  ) %>%
  mutate(across(starts_with("mean_"), ~round(. * 1e9, 4))) %>%
  mutate(across(starts_with("se_"),   ~round(. * 1e9, 4))) %>%
  select(Strain, mut_class,
         nonmono_mean = mean_rate_nonmono, nonmono_SEM = se_rate_nonmono,
         mono_mean    = mean_rate_mono,    mono_SEM    = se_rate_mono,
         total_mean   = mean_rate_total,   total_SEM   = se_rate_total)

cat("\nRate table (x10^9, SR, primary results):\n")
print(rate_table)
write.csv(rate_table, "rate_table_rajaei_plus2_olddenom_monoATGC_unitAssoc.csv", row.names=FALSE)
cat("Saved: rate_table_rajaei_plus2_olddenom_monoATGC_unitAssoc.csv\n")

# =============================================================================
# STEP 7: RATE RATIO with BOOTSTRAP 95% CI
# mu_mono / mu_nonmono for each class x strain
# This is the key Rajaei finding: AT->TA ratio ~7x, +/-1bp indels ~26-39x
#
# Bootstrap: resample lines with replacement within strain,
#            compute mean rate in each context, take ratio
# =============================================================================
cat("\nBootstrapping rate ratios...\n")

bootstrap_ratio <- function(rates_df, strain_f, class_f) {
  mono    <- rates_df %>%
    filter(Strain==strain_f, mut_class==class_f, context=="mono") %>%
    pull(rate)
  nonmono <- rates_df %>%
    filter(Strain==strain_f, mut_class==class_f, context=="nonmono") %>%
    pull(rate)

  n <- length(mono)
  if (n == 0 || mean(nonmono, na.rm=TRUE) == 0) return(NULL)

  # Point estimate
  obs_ratio <- mean(mono, na.rm=TRUE) / mean(nonmono, na.rm=TRUE)

  # Bootstrap
  boot_ratios <- numeric(N_BOOT)
  for (b in seq_len(N_BOOT)) {
    idx <- sample.int(n, n, replace=TRUE)
    m   <- mean(mono[idx],    na.rm=TRUE)
    nm  <- mean(nonmono[idx], na.rm=TRUE)
    boot_ratios[b] <- if (nm > 0) m / nm else NA_real_
  }
  ci <- quantile(boot_ratios, c(0.025, 0.975), na.rm=TRUE)

  data.frame(
    Strain       = strain_f,
    mut_class    = class_f,
    obs_ratio    = obs_ratio,
    ci_lo        = ci[1],
    ci_hi        = ci[2],
    log2_ratio   = log2(obs_ratio),
    n_lines      = n
  )
}

ratio_results <- bind_rows(lapply(c("HK104","PB800"), function(str) {
  bind_rows(lapply(ALL_CLASSES, function(cl) {
    bootstrap_ratio(rates_sr, str, cl)
  }))
})) %>%
  mutate(
    label = case_when(
      mut_class %in% SNP_CLASSES  ~
        SNP_LABELS[match(mut_class, SNP_CLASSES)],
      mut_class %in% INDEL_BINS   ~
        INDEL_LABELS[match(mut_class, INDEL_BINS)],
      TRUE ~ mut_class
    )
  )

write.csv(ratio_results, "rate_ratio_bootstrap_plus2_olddenom_monoATGC_unitAssoc.csv", row.names=FALSE)
cat("Saved: rate_ratio_bootstrap_plus2_olddenom_monoATGC_unitAssoc.csv\n")

cat("\nRate ratios (mu_mono / mu_nonmono, SR):\n")
print(ratio_results %>%
        filter(mut_class %in% c(SNP_CLASSES,"del_1","ins_1")) %>%
        select(Strain, label, obs_ratio, ci_lo, ci_hi) %>%
        arrange(Strain, obs_ratio))

# =============================================================================
# Define snv_clean here -- used by both the MC-FET (Step 8) and the
# spectrum figure (Step 9). Filters out any mutations where mono_repeat
# or Strain is NA (e.g. mutations on scaffolds absent from the repeat BED).
snv_clean <- sr %>%
  filter(TYPE=="snp",
         !is.na(mut_class),  mut_class %in% SNP_CLASSES,
         !is.na(mono_repeat), !is.na(Strain))

# =============================================================================
# STEP 8: MONTE CARLO FISHER'S EXACT TEST ON 6-CLASS SPECTRUM
# Rajaei et al. used MC-FET (10^6 iterations) to compare spectra
#
# Comparisons:
#   A. Mono vs nonmono spectrum -- within each strain (key context comparison)
#   B. HK104 vs PB800 -- within each context
# =============================================================================
cat("\nRunning Monte Carlo Fisher's Exact Tests...\n")

# Get SNV counts by strain x context (summed over lines)
# Use snv_clean (NA mono_repeat and NA Strain already removed)
# Same mutations as the spectrum figure -- ensures MC-FET tests the same data
snv_counts <- snv_clean %>%
  mutate(context = ifelse(mono_repeat, "mono", "nonmono")) %>%
  count(Strain, context, mut_class) %>%
  complete(Strain    = c("HK104","PB800"),
           context   = c("mono","nonmono"),
           mut_class = SNP_CLASSES,
           fill      = list(n=0))

mc_fet <- function(mat, B=1e6) {
  tryCatch(
    fisher.test(mat, simulate.p.value=TRUE, B=B)$p.value,
    error=function(e) NA_real_
  )
}

spectrum_tests <- bind_rows(
  # A. Mono vs nonmono within each strain
  bind_rows(lapply(c("HK104","PB800"), function(str) {
    mat <- snv_counts %>%
      filter(Strain==str) %>%
      pivot_wider(names_from=context, values_from=n) %>%
      select(mono, nonmono) %>%
      as.matrix()
    rownames(mat) <- SNP_CLASSES
    data.frame(
      comparison = "mono_vs_nonmono",
      Strain1    = str, Strain2 = str,
      context    = "both",
      p_value    = mc_fet(mat)
    )
  })),
  # B. HK104 vs PB800 within each context
  bind_rows(lapply(c("mono","nonmono","total"), function(ctx) {
    if (ctx=="total") {
      sub <- sr %>%
        filter(TYPE=="snp", !is.na(mut_class), mut_class %in% SNP_CLASSES) %>%
        count(Strain, mut_class) %>%
        complete(Strain, mut_class=SNP_CLASSES, fill=list(n=0))
    } else {
      sub <- snv_counts %>% filter(context==ctx)
    }
    mat <- sub %>%
      pivot_wider(names_from=Strain, values_from=n) %>%
      select(HK104, PB800) %>%
      as.matrix()
    rownames(mat) <- SNP_CLASSES
    data.frame(
      comparison = "HK104_vs_PB800",
      Strain1    = "HK104", Strain2="PB800",
      context    = ctx,
      p_value    = mc_fet(mat)
    )
  }))
)

cat("\nSpectrum comparison p-values (MC-FET, 10^6 iterations):\n")
print(spectrum_tests)
write.csv(spectrum_tests, "spectrum_stats_plus2_olddenom_monoATGC_unitAssoc.csv", row.names=FALSE)
cat("Saved: spectrum_stats_plus2_olddenom_monoATGC_unitAssoc.csv\n")

# =============================================================================
# STEP 9: FIGURES
# =============================================================================
cat("\nGenerating figures...\n")

# -- Fig 1: SNV rates by context (Rajaei Fig 2A-C equivalent) ---------------
# Build a single combined data frame and use facet_wrap(scales="fixed")
# so all three panels share the SAME y-axis -- this is how Rajaei Fig 2 works
# and is the only way to visually compare nonmono vs mono vs total rates.
# Y-axis is numerical (rates already scaled x10^9), not scientific notation.

make_rate_figure <- function(classes, labels, add_vline=FALSE) {
  d <- all_rate_summary %>%
    filter(dataset=="SR", mut_class %in% classes) %>%
    mutate(
      label   = factor(labels[match(mut_class, classes)], levels=labels),
      context = factor(context,
                       levels = c("nonmono","mono","total"),
                       labels = c("Non-mononucleotide","Mononucleotide","Total")),
      mean_r  = mean_rate * 1e9,
      se_r    = se_rate   * 1e9
    )

  p <- ggplot(d, aes(x=label, y=mean_r, fill=Strain)) +
    geom_col(position=position_dodge(0.7), width=0.65,
             color="white", linewidth=0.3) +
    geom_errorbar(
      aes(ymin=mean_r - se_r, ymax=mean_r + se_r),
      position=position_dodge(0.7), width=0.22, linewidth=0.5
    ) +
    scale_fill_manual(values=c(HK104=COL_HK, PB800=COL_PB)) +
    scale_y_continuous(
      # Numerical labels (not scientific) -- values are already x10^9
      labels = scales::number_format(accuracy=0.01),
      expand = expansion(mult=c(0, 0.15))
    ) +
    labs(x=NULL,
         y=expression(mu~(x10^{-9}~"per site per gen"))) +
    # scales="fixed" enforces the same y-axis range across all three panels
    facet_wrap(~context, nrow=1, scales="fixed") +
    theme_rajaei()

  if (add_vline) {
    p <- p + geom_vline(xintercept=5.5, linetype="dashed",
                        color="grey50", linewidth=0.5)
  }
  p
}

fig_snv <- make_rate_figure(SNP_CLASSES, SNP_LABELS) +
  plot_layout(guides="collect") & theme(legend.position="top", legend.title=element_text(size=18), legend.text=element_text(size=18))

save_fig(fig_snv, "Fig_SNV_rates_plus2_olddenom_monoATGC_unitAssoc", 14, 5)

# -- Fig 2: Indel rates by context (Rajaei Fig 2D-F equivalent) -------------
fig_ind <- make_rate_figure(INDEL_BINS, INDEL_LABELS, add_vline=TRUE) +
  plot_layout(guides="collect") & theme(legend.position="top", legend.title=element_text(size=18), legend.text=element_text(size=18))

save_fig(fig_ind, "Fig_Indel_rates_plus2_olddenom_monoATGC_unitAssoc", 16, 5)

# -- Fig 3: SNV SPECTRUM (proportions) by context (Rajaei Fig 4A-C equiv) ---
# Proportions: count of class / total SNVs in that context x strain
# This is NOT rates -- just frequency distribution (no normalisation by callable)
#
# Filter: drop any row where mono_repeat or Strain is NA before counting.
# The previous version had a case_when(is.na(mono_repeat) ~ "total") branch
# which caused NA-Strain rows to appear as a grey bar in the figure.

snv_spectrum <- bind_rows(
  # mono and nonmono contexts: use mono_repeat flag directly
  snv_clean %>% mutate(context = ifelse(mono_repeat, "mono", "nonmono")),
  # total context: all mutations regardless of repeat context
  snv_clean %>% mutate(context = "total")
) %>%
  count(Strain, context, mut_class) %>%
  group_by(Strain, context) %>%
  mutate(prop  = n / sum(n),
         label = factor(SNP_LABELS[match(mut_class, SNP_CLASSES)],
                        levels = SNP_LABELS)) %>%
  ungroup()

spec_panel <- function(ctx, subtitle_txt) {
  snv_spectrum %>%
    filter(context==ctx) %>%
    ggplot(aes(x=label, y=prop, fill=Strain)) +
      geom_col(position=position_dodge(0.7), width=0.65,
               color="white", linewidth=0.3) +
      scale_fill_manual(values=c(HK104=COL_HK, PB800=COL_PB)) +
      scale_y_continuous(labels=scales::percent,
                         expand=expansion(mult=c(0,0.12)),
                         limits=c(0,0.75)) +
      labs(x=NULL, y="Proportion", subtitle=subtitle_txt) +
      theme_rajaei()
}

fig_spec <- (spec_panel("nonmono","Non-mononucleotide") |
             spec_panel("mono",   "Mononucleotide")     |
             spec_panel("total",  "Total")) +
  plot_layout(guides="collect") & theme(legend.position="top", legend.title=element_text(size=18), legend.text=element_text(size=18))

save_fig(fig_spec, "Fig_SNV_spectrum_plus2_olddenom_monoATGC_unitAssoc", 14, 5)

# -- Fig 4: RATE RATIO with bootstrap 95% CI ---------------------------------
# mu_mono / mu_nonmono for each class; dashed line at 1 (no enrichment)
# This is the key summary figure -- Rajaei's ~7x AT->TA enrichment equivalent

ratio_plot_data <- ratio_results %>%
  mutate(
    class_type = ifelse(mut_class %in% SNP_CLASSES, "SNV", "Indel"),
    label = factor(label,
                   levels=c(SNP_LABELS, INDEL_LABELS))
  )

fig_ratio <- ratio_plot_data %>%
  ggplot(aes(x=label, y=obs_ratio,
             ymin=ci_lo, ymax=ci_hi, color=Strain)) +
    geom_hline(yintercept=1, linetype="dashed",
               color="grey40", linewidth=0.6) +
    geom_errorbar(position=position_dodge(0.5),
                  width=0.3, linewidth=0.7) +
    geom_point(position=position_dodge(0.5),
               size=2.5, shape=19) +
    scale_color_manual(values=c(HK104=COL_HK, PB800=COL_PB)) +
    scale_y_log10(
      breaks=c(0.1,0.5,1,2,5,10,20,50),
      labels=c("0.1","0.5","1","2","5","10","20","50")
    ) +
    labs(x=NULL,
         y=expression(mu["mono"]/mu["nonmono"]~"(log"[10]~"scale)")) +
    facet_wrap(~class_type, scales="free_x") +
    theme_rajaei() +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=18))

save_fig(fig_ratio, "Fig_rate_ratio_plus2_olddenom_monoATGC_unitAssoc", 13, 5)

# -- Fig 5: Ts/Tv IN REPEAT vs NON-REPEAT ------------------------------------
# Filter: remove NA Strain (mutations on scaffolds absent from repeat BED)
#         and NA mono_repeat (same cause) before counting
tstv_data <- sr %>%
  filter(TYPE=="snp",
         !is.na(mut_class),    mut_class %in% SNP_CLASSES,
         !is.na(mono_repeat),
         !is.na(Strain),       Strain %in% c("HK104","PB800")) %>%
  mutate(
    context = ifelse(mono_repeat, "Mononucleotide", "Non-mononucleotide"),
    context = factor(context, levels=c("Mononucleotide","Non-mononucleotide")),
    ts_tv   = ifelse(mut_class %in% TS_CLASSES, "Ts", "Tv")
  ) %>%
  count(Strain, context, ts_tv) %>%
  pivot_wider(names_from=ts_tv, values_from=n, values_fill=0) %>%
  mutate(TsTv = Ts / Tv)

cat("\nTs/Tv by context:\n")
print(tstv_data)

fig_tstv <- tstv_data %>%
  ggplot(aes(x=context, y=TsTv, fill=Strain)) +
    geom_col(position=position_dodge(0.7), width=0.6,
             color="white", linewidth=0.3) +
    scale_fill_manual(values=c(HK104=COL_HK, PB800=COL_PB),
                      limits=c("HK104","PB800")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
    labs(x=NULL, y="Ts/Tv ratio") +
    theme_rajaei()

save_fig(fig_tstv, "Fig_TsTv_context_plus2_olddenom_monoATGC_unitAssoc", 6, 4)

# =============================================================================
# EXTENSION: Per-unit enrichment heatmap (beyond Rajaei) -- corrected null
# For each SNV class x repeat unit (A/T/C/G), the null is conditioned on the
# mutable base context of that SNV class:
#   GC_to_* classes are tested against C/G-repeat bp among total GC bp.
#   AT_to_* classes are tested against A/T-repeat bp among total AT bp.
# This avoids testing AT-context mutations against the whole genome, which would
# be an invalid null for base-substitution classes.
# =============================================================================
cat("\n--- Extension: Per-unit enrichment (context-corrected null; beyond Rajaei) ---\n")

unit_bp_vec <- setNames(unit_bp$unit_bp, unit_bp$repeat_unit)
unit_bp_vec[setdiff(UNIT_LEVELS, names(unit_bp_vec))] <- 0
unit_bp_vec <- unit_bp_vec[UNIT_LEVELS]

unit_null_for_class <- function(class_f, unit_f) {
  # Heatmap now means: distribution across repeat-associated units
  # among mutations inside the plus2 mono context. Therefore all units
  # are admissible, and the null is the A/T/C/G composition of the
  # repeat-associated intervals, not source-base compatibility.
  val <- unit_bp_assoc$p_null_assoc[unit_bp_assoc$repeat_unit == unit_f]
  if (length(val) == 0) return(NA_real_)
  as.numeric(val[1])
}

unit_enr <- bind_rows(lapply(c("HK104", "PB800"), function(str) {
  snv <- sr %>%
    filter(Strain == str,
           TYPE == "snp",
           mono_repeat == TRUE,
           !is.na(mut_class),
           mut_class %in% SNP_CLASSES,
           !is.na(repeat_unit))

  bind_rows(lapply(SNP_CLASSES, function(cl) {
    sub   <- snv %>% filter(mut_class == cl)
    n_tot <- nrow(sub)

    allowed_units <- UNIT_LEVELS

    bind_rows(lapply(allowed_units, function(unit) {
      p_u  <- unit_null_for_class(cl, unit)
      if (is.na(p_u) || p_u <= 0 || n_tot == 0) return(NULL)
      n_in <- sum(sub$repeat_unit == unit, na.rm = TRUE)
      bt   <- binom.test(n_in, n_tot, p = p_u, alternative = "two.sided")
      data.frame(
        Strain = str, class = cl,
        label = SNP_LABELS[match(cl, SNP_CLASSES)],
        unit = unit, n_total = n_tot, n_in_unit = n_in,
        obs_frac = n_in / n_tot, null_frac = p_u,
        enrichment_ratio = (n_in / n_tot) / p_u,
        p_value = bt$p.value, p_adj = NA, significant = FALSE
      )
    }))
  }))
})) %>%
  group_by(Strain) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         significant = !is.na(p_adj) & p_adj < 0.05) %>%
  ungroup()

write.csv(unit_enr, "unit_enrichment_extension_plus2_olddenom_monoATGC_unitAssoc.csv", row.names = FALSE)

fig_heatmap <- unit_enr %>%
  mutate(
    label    = factor(label, levels = SNP_LABELS),
    unit     = factor(unit,  levels = UNIT_LEVELS),
    log2_enr = log2(enrichment_ratio),
    sig_lbl  = ifelse(significant, "*", "")
  ) %>%
  ggplot(aes(x = unit, y = label, fill = log2_enr)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_lbl), size = 5, vjust = 0.8) +
    scale_fill_gradient2(
      low = "#1F77B4", mid = "white", high = "#D62728", midpoint = 0,
      name = "log2(enrichment)\n* FDR<0.05"
    ) +
    facet_wrap(~Strain, nrow = 1) +
    labs(x = "Repeat-associated unit (tract +/- 2 bp)", y = "SNV class") +
    theme_rajaei() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=18),
          panel.grid = element_blank())

save_fig(fig_heatmap, "Fig_unit_heatmap_plus2_olddenom_monoATGC_unitAssoc", 9, 5)


# =============================================================================
# EXTENSION v2: Improved per-unit enrichment heatmap (Option A + selective B)
# Adds:
#   - n=X count annotation inside every tested cell (so viewer sees count base)
#   - distinct visual fills for three cell categories:
#       * tested  -> color by log2(enrichment)
#       * n=0     -> medium grey
#       * impossible -> very light grey + em-dash
#   - thick horizontal divider between AT-class (top 3) and GC-class (bottom 3)
#     rows to make the structural split obvious
#   - dotted border overlay on cells with n < 10 (visual "low count" warning)
#   - prominent bold asterisks for FDR<0.05 cells
# Keeps the original Fig_unit_heatmap_plus2_olddenom_monoATGC_unitAssoc.* output as the legacy version;
# this v2 produces Fig_unit_heatmap_plus2_olddenom_monoATGC_unitAssoc_v2.*
# =============================================================================
cat("\n--- Extension v2: Improved per-unit heatmap (count-annotated, count-binned) ---\n")

# 1. Build complete grid including mechanically-impossible combinations
all_combos <- expand.grid(
  Strain = c("HK104", "PB800"),
  class  = SNP_CLASSES,
  unit   = UNIT_LEVELS,
  stringsAsFactors = FALSE
) %>%
  mutate(
    label      = SNP_LABELS[match(class, SNP_CLASSES)],
    base_class = ifelse(class %in% c("GC_to_AT", "GC_to_CG", "GC_to_TA"), "GC", "AT"),
    # With repeat-associated regions (+/-2 bp), every SNV class can occur
    # around any repeat unit because flanking bases need not match the run base.
    compatible = TRUE
  )

# 2. Merge with computed enrichments
plot_df_v2 <- all_combos %>%
  left_join(unit_enr %>% select(Strain, class, unit,
                                n_in_unit, enrichment_ratio, p_adj, significant),
            by = c("Strain", "class", "unit")) %>%
  mutate(
    cell_type = case_when(
      !compatible                       ~ "impossible",
      is.na(n_in_unit) | n_in_unit == 0 ~ "zero_obs",
      TRUE                              ~ "tested"
    ),
    log2_enr  = ifelse(cell_type == "tested" & enrichment_ratio > 0,
                       log2(enrichment_ratio), NA_real_),
    low_n     = !is.na(n_in_unit) & n_in_unit > 0 & n_in_unit < 10,
    sig_lbl   = ifelse(!is.na(significant) & significant, "*", ""),
    n_lbl     = case_when(
      cell_type == "impossible" ~ "\u2014",       # em-dash
      cell_type == "zero_obs"   ~ "n=0",
      TRUE                      ~ paste0("n=", n_in_unit)
    ),
    label     = factor(label, levels = rev(SNP_LABELS)),  # AT-classes appear at top
    unit      = factor(unit,  levels = UNIT_LEVELS)
  )

# 3. Build the heatmap (4 layered geom_tile calls for the 3 cell categories
#    plus a low-n border overlay)
fig_heatmap_v2 <- ggplot(plot_df_v2, aes(x = unit, y = label)) +
  # Layer 1: impossible cells (very light grey, em-dash inside)
  geom_tile(data = function(d) filter(d, cell_type == "impossible"),
            fill = "#f0f0f0", color = "#dddddd", linewidth = 0.4) +
  # Layer 2: zero-obs cells (medium grey)
  geom_tile(data = function(d) filter(d, cell_type == "zero_obs"),
            fill = "#bcbcbc", color = "white", linewidth = 1.0) +
  # Layer 3: tested cells (red/blue gradient by log2 enrichment)
  geom_tile(data = function(d) filter(d, cell_type == "tested"),
            aes(fill = log2_enr), color = "white", linewidth = 1.0) +
  # Layer 4: dotted border for low-count cells
  geom_tile(data = function(d) filter(d, low_n),
            fill = NA, color = "#444444", linewidth = 0.9, linetype = "dotted") +
  # Asterisks for FDR-significant cells
  geom_text(aes(label = sig_lbl), size = 7, vjust = -0.45, fontface = "bold") +
  # n=X counts (small, below center)
  geom_text(aes(label = n_lbl), size = 4.5, vjust = 1.4, color = "#1a1a1a") +
  # Horizontal divider between AT (top 3 rows) and GC (bottom 3 rows)
  geom_hline(yintercept = 3.5, color = "#333333", linewidth = 0.9) +
  scale_fill_gradient2(
    low = "#1F77B4", mid = "white", high = "#D62728", midpoint = 0,
    limits = c(-1.5, 1.5), oob = scales::squish,
    name = expression(log[2] ~ "(enrichment)")
  ) +
  facet_wrap(~ Strain, nrow = 1) +
  labs(x = "Repeat-associated unit (tract +/- 2 bp)", y = "SNV class",
       caption = paste("* = FDR < 0.05 (BH-adjusted)  |  grey = n=0 (tested, no events)  |  ",
                       "dotted border = low count (n < 10)")) +
  theme_rajaei() +
  theme(
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    panel.grid   = element_blank(),
    plot.caption = element_text(size = 13, hjust = 0, color = "#555555",
                                margin = margin(t = 8)),
    legend.position = "right",
    strip.text   = element_text(size = 18, face = "bold")
  )

save_fig(fig_heatmap_v2, "Fig_unit_heatmap_plus2_olddenom_monoATGC_unitAssoc_v2", 11.5, 5)

cat("Improved per-unit heatmap saved as Fig_unit_heatmap_plus2_olddenom_monoATGC_unitAssoc_v2.pdf/.png\n")


# =============================================================================
# DONE
# =============================================================================
cat("\n=== Complete. All outputs saved. ===\n")
cat("\n-- Rajaei-framework outputs: --\n")
cat("  rate_table_rajaei_plus2_olddenom_monoATGC_unitAssoc.csv      (Table 1 equivalent)\n")
cat("  rate_ratio_bootstrap_plus2_olddenom_monoATGC_unitAssoc.csv   (mu_mono/mu_nonmono + 95% CI)\n")
cat("  spectrum_stats_plus2_olddenom_monoATGC_unitAssoc.csv         (MC-FET p-values)\n")
cat("  Fig_SNV_rates_plus2_olddenom_monoATGC_unitAssoc.pdf/.png     (Fig 2A-C equivalent)\n")
cat("  Fig_Indel_rates_plus2_olddenom_monoATGC_unitAssoc.pdf/.png   (Fig 2D-F equivalent)\n")
cat("  Fig_SNV_spectrum_plus2_olddenom_monoATGC_unitAssoc.pdf/.png  (Fig 4A-C equivalent)\n")
cat("  Fig_rate_ratio_plus2_olddenom_monoATGC_unitAssoc.pdf/.png    (summary: ratio + CI)\n")
cat("  Fig_TsTv_context_plus2_olddenom_monoATGC_unitAssoc.pdf/.png\n")
cat("\n-- Extension output: --\n")
cat("  unit_enrichment_extension_plus2_olddenom_monoATGC_unitAssoc.csv\n")
cat("  Fig_unit_heatmap_plus2_olddenom_monoATGC_unitAssoc.pdf/.png      (A/T/C/G run context, original)\n")
cat("  Fig_unit_heatmap_plus2_olddenom_monoATGC_unitAssoc_v2.pdf/.png   (A/T/C/G run context, improved with count annotations)\n")

cat("\n-- Key numbers to compare with Rajaei Table 1 (C. elegans N2): --\n")
cat("  AT->TA nonmono: ~0.44 x10^-9 | mono: ~3.34 x10^-9 | ratio: ~7.6x\n")
cat("  Del-1  nonmono: ~0.19 x10^-9 | mono: ~1.43 x10^-9 | ratio: ~7.5x\n")
cat("  Ins-1  nonmono: ~0.07 x10^-9 | mono: ~1.07 x10^-9 | ratio: ~15x\n")
cat("  Your briggsae values are in rate_ratio_bootstrap.csv\n")

