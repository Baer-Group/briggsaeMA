# =============================================================================
# Rate Table -- C. briggsae MA Analysis (HK104 + PB800)
# =============================================================================
# Inputs:
#   master_table.csv          -- per-line mutation counts, callable sites,
#                                generations, Strain
#   recall_summary_grand.xlsx -- FTR_rate per sample per replicate
#
# Correction: FtR only (no FP subtraction -- FP handled upstream via
#             replicate concordance and IGV validation)
#
# FtR rate source: B1_SNP_DP3, B2_SNP_DP3, B1_Indel_DP3, B2_Indel_DP3
#   HK104 lines: sample number <= 298
#   PB800 lines: sample number >= 300
#   Pool B1 + B2 per strain; use pre-computed FTR_rate column directly
#
# Rate formula:
#   Step 1 -- corrected count: observed + Poisson(observed x p_ftr)
#   Step 2 -- rate:            corrected / (callable_sites x generations)
#   NOTE: callable_sites and FtR_rate are NEVER used in the same formula term
#
# Output:
#   rate_table_simple.csv     -- mu_bar (mean +/- SEM per strain)
#   rate_table_bootstrap.csv  -- mu_hat (mean, 95% CI after FtR correction)
#   rate_table_display.csv    -- combined display table matching Table 2-1 format
# =============================================================================

library(rstudioapi)
library(dplyr)
library(readxl)

setwd(dirname(getActiveDocumentContext()$path))
cat("Working directory:", getwd(), "\n")

# =============================================================================
# USER-ADJUSTABLE COLUMN NAMES
# Adjust these if your master_table column names differ
# =============================================================================
GENOME_SIZE  <- 106196309   # C. briggsae reference genome (bp)
N_BOOT       <- 10000
set.seed(2025)

# master_table column names
COL_STRAIN   <- "Strain"       # values: "HK104", "PB800"
COL_LINE     <- "Line"
COL_GEN      <- "Generations"  # per-line generation count
COL_CALLABLE <- "SR_callable_3x"  # callable sites at 3x threshold

# SNP count columns (6 classes -- summed to get total SNV)
SNP_COLS <- c("SR_3x_GC_to_AT", "SR_3x_GC_to_CG", "SR_3x_GC_to_TA",
              "SR_3x_AT_to_GC", "SR_3x_AT_to_CG", "SR_3x_AT_to_TA")

# Indel count columns -- deletions and insertions separately
DEL_COLS <- c("SR_3x_del_10plus", "SR_3x_del_6_10", "SR_3x_del_3_5",
              "SR_3x_del_2",      "SR_3x_del_1")
INS_COLS <- c("SR_3x_ins_1",      "SR_3x_ins_2",   "SR_3x_ins_3_5",
              "SR_3x_ins_6_10",   "SR_3x_ins_10plus")

# =============================================================================
# STEP 1: EXTRACT p_ftr FROM RECALL FILE
# Sheets used: B1_SNP_DP3, B2_SNP_DP3, B1_Indel_DP3, B2_Indel_DP3
# HK104: sample number <= 298
# PB800: sample number >= 300
# Pool B1 + B2 per strain; mean FTR_rate = p_ftr
# =============================================================================
cat("\n--- Extracting FtR rates from recall file ---\n")

recall_file <- "recall_summary_grand.xlsx"

read_recall_sheet <- function(sheet_name) {
  df <- read_excel(recall_file, sheet = sheet_name)
  # Remove summary rows (Mean, Median) -- keep only sample rows
  df <- df %>% filter(!is.na(SAMPLE), !SAMPLE %in% c("Mean", "Median"))
  # Extract numeric prefix from sample name (e.g. "202_POOLRET91..." -> 202)
  df$sample_num <- as.integer(sub("_.*", "", df$SAMPLE))
  # Assign strain by numeric cutoff
  df$strain <- ifelse(df$sample_num <= 298, "HK104", "PB800")
  # Keep only needed columns to avoid type conflicts across sheets
  # (Recall_rate and FN_rate can be double in one sheet, character in another)
  df %>%
    select(SAMPLE, sample_num, strain, FTR_rate) %>%
    mutate(FTR_rate = as.numeric(FTR_rate))
}

# Read the four DP3 sheets
b1_snp  <- read_recall_sheet("B1_SNP_DP3")
b2_snp  <- read_recall_sheet("B2_SNP_DP3")
b1_ind  <- read_recall_sheet("B1_Indel_DP3")
b2_ind  <- read_recall_sheet("B2_Indel_DP3")

# Pool B1 + B2 per mutation type, then get mean FTR_rate per strain
get_p_ftr <- function(df1, df2, type_label) {
  pooled <- bind_rows(df1, df2)
  result <- pooled %>%
    group_by(strain) %>%
    summarise(
      p_ftr      = mean(FTR_rate, na.rm = TRUE),
      n_samples  = n(),
      .groups    = "drop"
    ) %>%
    mutate(type = type_label)
  result
}

p_ftr_snp   <- get_p_ftr(b1_snp, b2_snp, "SNP")
p_ftr_indel <- get_p_ftr(b1_ind, b2_ind, "Indel")

cat("FtR rates extracted:\n")
print(bind_rows(p_ftr_snp, p_ftr_indel))

# Extract scalar p_ftr values for each strain x type
p_ftr <- list(
  HK104_SNP   = p_ftr_snp   %>% filter(strain == "HK104") %>% pull(p_ftr),
  PB800_SNP   = p_ftr_snp   %>% filter(strain == "PB800") %>% pull(p_ftr),
  HK104_Indel = p_ftr_indel %>% filter(strain == "HK104") %>% pull(p_ftr),
  PB800_Indel = p_ftr_indel %>% filter(strain == "PB800") %>% pull(p_ftr)
)
cat("\nFtR rate summary:\n")
for (nm in names(p_ftr)) cat(sprintf("  p_ftr[%s] = %.4f\n", nm, p_ftr[[nm]]))

# =============================================================================
# STEP 2: LOAD MASTER TABLE AND COMPUTE PER-LINE COUNTS
# =============================================================================
cat("\n--- Loading master table ---\n")

mt <- read.csv("master_table.csv", stringsAsFactors = FALSE)
cat(sprintf("  %d lines loaded\n", nrow(mt)))

# Compute aggregate counts per line
mt <- mt %>%
  filter(!is.na(.data[[COL_STRAIN]]),
         !is.na(.data[[COL_CALLABLE]]),
         !is.na(.data[[COL_GEN]])) %>%
  mutate(
    count_SNV   = rowSums(across(all_of(SNP_COLS)), na.rm = TRUE),
    count_DEL   = rowSums(across(all_of(DEL_COLS)), na.rm = TRUE),
    count_INS   = rowSums(across(all_of(INS_COLS)), na.rm = TRUE),
    count_TOTAL = count_SNV + count_DEL + count_INS,
    callable    = .data[[COL_CALLABLE]],
    gens        = .data[[COL_GEN]],
    strain      = .data[[COL_STRAIN]]
  ) %>%
  select(Line = all_of(COL_LINE), strain, callable, gens,
         count_SNV, count_DEL, count_INS, count_TOTAL)

cat(sprintf("  HK104: %d lines | PB800: %d lines\n",
            sum(mt$strain == "HK104"), sum(mt$strain == "PB800")))

# =============================================================================
# HELPERS
# =============================================================================
sem    <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
ci95   <- function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)

# Denominator: callable_sites x generations (never mixed with FtR)
denom  <- function(callable, gens) as.numeric(callable) * as.numeric(gens)

# FtR-corrected count: observed + Poisson(observed x p_ftr)
# p_ftr and counts are used SEPARATELY -- p_ftr never enters the denominator
ftr_correct <- function(counts, p_ftr_val) {
  counts + rpois(length(counts), lambda = counts * p_ftr_val)
}

# Per-line rate: corrected_count / (callable x generations)
line_rate <- function(counts, callable, gens) counts / denom(callable, gens)

# =============================================================================
# STEP 3: OBSERVED RATES (mu_bar) -- no correction
# =============================================================================
cat("\n--- Computing observed rates (mu_bar) ---\n")

mut_classes   <- c("SNV", "DEL", "INS", "TOTAL")
count_col_map <- c(SNV   = "count_SNV",
                   DEL   = "count_DEL",
                   INS   = "count_INS",
                   TOTAL = "count_TOTAL")

simple_out <- bind_rows(lapply(c("HK104", "PB800"), function(str) {
  df <- mt %>% filter(strain == str)
  bind_rows(lapply(mut_classes, function(cl) {
    r <- line_rate(df[[count_col_map[[cl]]]], df$callable, df$gens)
    data.frame(
      strain     = str,
      class      = cl,
      N_lines    = nrow(df),
      pct_genome = mean(df$callable / GENOME_SIZE * 100, na.rm = TRUE),
      pct_genome_sem = sem(df$callable / GENOME_SIZE * 100),
      mean_gens  = mean(df$gens, na.rm = TRUE),
      gens_sem   = sem(df$gens),
      mu_bar     = mean(r, na.rm = TRUE),
      mu_bar_sem = sem(r)
    )
  }))
}))

print(simple_out)
write.csv(simple_out, "rate_table_simple.csv", row.names = FALSE)
cat("Saved: rate_table_simple.csv\n")

# =============================================================================
# STEP 4: BOOTSTRAP WITH FtR CORRECTION (mu_hat)
# For each bootstrap replicate:
#   1. Resample lines with replacement (within strain)
#   2. Apply FtR correction to counts (step 1: numerator only)
#   3. Compute rate using callable x generations (step 2: denominator only)
# =============================================================================
cat("\n--- Running bootstrap (N =", N_BOOT, ") ---\n")

run_bootstrap <- function(df, strain_name) {
  n  <- nrow(df)
  p_snp   <- p_ftr[[paste0(strain_name, "_SNP")]]
  p_indel <- p_ftr[[paste0(strain_name, "_Indel")]]

  # Storage: rows = bootstrap replicates, cols = mutation classes
  boot_mat <- matrix(NA_real_, nrow = N_BOOT, ncol = 4,
                     dimnames = list(NULL, mut_classes))

  for (b in seq_len(N_BOOT)) {
    idx  <- sample.int(n, size = n, replace = TRUE)
    samp <- df[idx, ]

    # Step 1: FtR correction applied to counts (numerator)
    # NOTE: p_ftr is NEVER applied to callable or gens (the denominator)
    c_snv   <- ftr_correct(samp$count_SNV,   p_snp)
    c_del   <- ftr_correct(samp$count_DEL,   p_indel)
    c_ins   <- ftr_correct(samp$count_INS,   p_indel)
    c_total <- c_snv + c_del + c_ins

    # Step 2: rate = corrected_count / (callable x generations)
    den <- denom(samp$callable, samp$gens)
    boot_mat[b, "SNV"]   <- mean(c_snv   / den, na.rm = TRUE)
    boot_mat[b, "DEL"]   <- mean(c_del   / den, na.rm = TRUE)
    boot_mat[b, "INS"]   <- mean(c_ins   / den, na.rm = TRUE)
    boot_mat[b, "TOTAL"] <- mean(c_total / den, na.rm = TRUE)
  }
  boot_mat
}

# Run per strain
boot_HK <- run_bootstrap(mt %>% filter(strain == "HK104"), "HK104")
boot_PB <- run_bootstrap(mt %>% filter(strain == "PB800"), "PB800")
cat("Bootstrap complete.\n")

# Summarise bootstrap results
boot_summary <- bind_rows(lapply(list(HK104 = boot_HK, PB800 = boot_PB),
                                 function(bmat) {
  bind_rows(lapply(mut_classes, function(cl) {
    x  <- bmat[, cl]
    ci <- ci95(x)
    data.frame(class  = cl,
               mu_hat = mean(x, na.rm = TRUE),
               ci_lo  = ci[1],
               ci_hi  = ci[2])
  }))
}), .id = "strain")

write.csv(boot_summary, "rate_table_bootstrap.csv", row.names = FALSE)
cat("Saved: rate_table_bootstrap.csv\n")

# =============================================================================
# STEP 5: COMBINED DISPLAY TABLE (matches Table 2-1 format)
# Columns: Strain | N_Lines | %Genome (SEM) | t-bar (SEM) |
#          mu_SNV (SEM) | mu_DEL (SEM) | mu_INS (SEM) | mu_TOTAL (SEM) |
#          mu_hat_SNV (95% CI) | mu_hat_DEL (95% CI) |
#          mu_hat_INS (95% CI) | mu_hat_TOTAL (95% CI)
# Rates reported x 10^9
# =============================================================================
cat("\n--- Building display table ---\n")

scale <- 1e9  # report as x 10^9

display <- bind_rows(lapply(c("HK104", "PB800"), function(str) {
  s_row <- simple_out  %>% filter(strain == str)
  b_row <- boot_summary %>% filter(strain == str)

  # Metadata (same across all classes)
  meta <- s_row %>% filter(class == "SNV") %>%
    select(N_lines, pct_genome, pct_genome_sem, mean_gens, gens_sem)

  bind_rows(lapply(mut_classes, function(cl) {
    s <- s_row  %>% filter(class == cl)
    b <- b_row  %>% filter(class == cl)
    data.frame(
      Strain          = str,
      Class           = cl,
      N_Lines         = meta$N_lines,
      pct_Genome      = round(meta$pct_genome,     3),
      pct_Genome_SEM  = round(meta$pct_genome_sem, 3),
      t_bar           = round(meta$mean_gens,       1),
      t_bar_SEM       = round(meta$gens_sem,        2),
      mu_bar_x1e9     = round(s$mu_bar     * scale, 4),
      mu_bar_SEM_x1e9 = round(s$mu_bar_sem * scale, 4),
      mu_hat_x1e9     = round(b$mu_hat * scale, 4),
      mu_hat_CI_lo    = round(b$ci_lo  * scale, 4),
      mu_hat_CI_hi    = round(b$ci_hi  * scale, 4)
    )
  }))
}))

print(display, digits = 4)
write.csv(display, "rate_table_display.csv", row.names = FALSE)
cat("Saved: rate_table_display.csv\n")

cat("\n=== Done. Files written: ===\n")
cat("  rate_table_simple.csv\n")
cat("  rate_table_bootstrap.csv\n")
cat("  rate_table_display.csv\n")
