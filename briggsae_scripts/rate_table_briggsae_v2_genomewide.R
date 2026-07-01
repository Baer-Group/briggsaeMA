# =============================================================================
# Rate Table -- C. briggsae MA Analysis (HK104 + PB800)
# v2 -- GENOME-WIDE FtR MODEL
# =============================================================================
# Updated 2026-05-24 based on clarification that FtR is GENOME-WIDE failure:
#   Dummy mutations were inserted across the entire reference genome (proportional
#   to chromosome length). FtR = 1 - (recovered / total inserted). This captures
#   losses from BOTH (a) mutations falling in uncallable regions (no coverage) and
#   (b) mutations in callable regions but called as homozygous reference.
#
# Generative model under genome-wide FtR:
#     m_obs  ~  Poisson(mu × L_total × t × (1 - p_ftr))
#
# Implied per-bp mutation rate estimator:
#     mu_hat = m_obs / [L_total × t × (1 - p_ftr)]
#
# CHANGES FROM v1:
#   1. DENOMINATOR:
#         v1: L_callable (line-specific callable size at DP >= 3)
#         v2: L_total    (reference genome length, constant = 106,196,309 bp)
#      The corrected rate is now per bp of REFERENCE GENOME -- the biologically
#      meaningful per-bp mutation rate.
#
#   2. FtR CORRECTION:
#         v1: m_hat = m + Pois(m × p_ftr)        [E ~ m(1+p), additive, first-order]
#         v2: m_hat = m + Pois(m × p_ftr/(1-p))  [E   = m/(1-p), exact multiplicative
#                                                 expectation, with Poisson noise]
#      The exact inverse-recall form is mathematically correct under the genome-
#      wide model. Setting lambda = m × p/(1-p) makes E[m_hat] = m / (1 - p) exact.
#
#   3. OUTPUT:
#      For comparison and back-checking, the simple table reports BOTH:
#         mu_bar_per_callable (old denominator) -- not corrected
#         mu_bar_per_total    (new denominator) -- not corrected
#      The bootstrap output reports mu_hat in the new (per-reference-bp) convention.
#
# Inputs:
#   master_table.csv          -- per-line mutation counts, callable sites,
#                                generations, Strain
#   recall_summary_grand.xlsx -- FTR_rate per sample per replicate
#
# Outputs:
#   rate_table_simple.csv     -- mu_bar (mean +/- SEM per strain, both conventions)
#   rate_table_bootstrap.csv  -- mu_hat (mean, 95% CI; per reference bp)
#   rate_table_display.csv    -- combined display table
# =============================================================================

library(rstudioapi)
library(dplyr)
library(readxl)

setwd(dirname(getActiveDocumentContext()$path))
cat("Working directory:", getwd(), "\n")

# =============================================================================
# USER-ADJUSTABLE CONSTANTS
# =============================================================================
GENOME_SIZE  <- 106196309   # C. briggsae reference genome (bp) -- L_total
N_BOOT       <- 10000
set.seed(2025)

# master_table column names
COL_STRAIN   <- "Strain"        # values: "HK104", "PB800"
COL_LINE     <- "Line"
COL_GEN      <- "Generations"   # per-line generation count
COL_CALLABLE <- "SR_callable_3x"  # per-line callable sites at DP >= 3

# SNP count columns (6 classes -- summed to get total SNV)
SNP_COLS <- c("SR_3x_GC_to_AT", "SR_3x_GC_to_CG", "SR_3x_GC_to_TA",
              "SR_3x_AT_to_GC", "SR_3x_AT_to_CG", "SR_3x_AT_to_TA")

# Indel count columns
DEL_COLS <- c("SR_3x_del_10plus", "SR_3x_del_6_10", "SR_3x_del_3_5",
              "SR_3x_del_2",      "SR_3x_del_1")
INS_COLS <- c("SR_3x_ins_1",      "SR_3x_ins_2",   "SR_3x_ins_3_5",
              "SR_3x_ins_6_10",   "SR_3x_ins_10plus")

# =============================================================================
# STEP 1: EXTRACT p_ftr FROM RECALL FILE (genome-wide failure rate)
# Sheets used: B1_SNP_DP3, B2_SNP_DP3, B1_Indel_DP3, B2_Indel_DP3
# HK104: sample number <= 298
# PB800: sample number >= 300
# Pool B1 + B2 per strain; mean FTR_rate = p_ftr (genome-wide)
# =============================================================================
cat("\n--- Extracting FtR rates from recall file ---\n")
cat("    NOTE: FtR is GENOME-WIDE (dummies inserted across entire reference)\n")

recall_file <- "recall_summary_grand.xlsx"

read_recall_sheet <- function(sheet_name) {
  df <- read_excel(recall_file, sheet = sheet_name)
  df <- df %>% filter(!is.na(SAMPLE), !SAMPLE %in% c("Mean", "Median"))
  df$sample_num <- as.integer(sub("_.*", "", df$SAMPLE))
  df$strain <- ifelse(df$sample_num <= 298, "HK104", "PB800")
  df %>%
    select(SAMPLE, sample_num, strain, FTR_rate) %>%
    mutate(FTR_rate = as.numeric(FTR_rate))
}

b1_snp  <- read_recall_sheet("B1_SNP_DP3")
b2_snp  <- read_recall_sheet("B2_SNP_DP3")
b1_ind  <- read_recall_sheet("B1_Indel_DP3")
b2_ind  <- read_recall_sheet("B2_Indel_DP3")

get_p_ftr <- function(df1, df2, type_label) {
  pooled <- bind_rows(df1, df2)
  pooled %>%
    group_by(strain) %>%
    summarise(
      p_ftr      = mean(FTR_rate, na.rm = TRUE),
      p_ftr_sd   = sd(FTR_rate, na.rm = TRUE),
      n_samples  = n(),
      .groups    = "drop"
    ) %>%
    mutate(type = type_label)
}

p_ftr_snp   <- get_p_ftr(b1_snp, b2_snp, "SNP")
p_ftr_indel <- get_p_ftr(b1_ind, b2_ind, "Indel")

cat("FtR rates extracted (genome-wide):\n")
print(bind_rows(p_ftr_snp, p_ftr_indel))

# Extract scalar p_ftr values
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
# HELPERS  (UPDATED FOR GENOME-WIDE FtR)
# =============================================================================
sem  <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
ci95 <- function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)

# CORRECTED FtR correction (genome-wide): m_hat = m + Pois(m × p/(1-p))
#   E[m_hat]      = m / (1 - p)              [exact multiplicative]
#   Var[m_hat | m, p] adds Poisson sim noise as in v1
# To switch to a deterministic correction, set m_hat = m / (1 - p) instead.
ftr_correct <- function(counts, p_ftr_val) {
  lambda <- counts * p_ftr_val / (1 - p_ftr_val)
  counts + rpois(length(counts), lambda = lambda)
}

# Denominators (no per-line FtR, no callable-size correction in the new convention)
denom_callable <- function(callable, gens) as.numeric(callable) * as.numeric(gens)
denom_total    <- function(gens)           GENOME_SIZE * as.numeric(gens)

# =============================================================================
# STEP 3: OBSERVED RATES (mu_bar) -- both conventions, no FtR correction
# =============================================================================
cat("\n--- Computing observed rates (mu_bar) -- both conventions ---\n")

mut_classes   <- c("SNV", "DEL", "INS", "TOTAL")
count_col_map <- c(SNV   = "count_SNV",
                   DEL   = "count_DEL",
                   INS   = "count_INS",
                   TOTAL = "count_TOTAL")

simple_out <- bind_rows(lapply(c("HK104", "PB800"), function(str) {
  df <- mt %>% filter(strain == str)
  bind_rows(lapply(mut_classes, function(cl) {
    counts   <- df[[count_col_map[[cl]]]]
    r_callab <- counts / denom_callable(df$callable, df$gens)
    r_total  <- counts / denom_total(df$gens)
    data.frame(
      strain                  = str,
      class                   = cl,
      N_lines                 = nrow(df),
      pct_genome              = mean(df$callable / GENOME_SIZE * 100, na.rm = TRUE),
      pct_genome_sem          = sem(df$callable / GENOME_SIZE * 100),
      mean_gens               = mean(df$gens, na.rm = TRUE),
      gens_sem                = sem(df$gens),
      # Old convention (per callable bp): kept for back-compatibility
      mu_bar_per_callable     = mean(r_callab, na.rm = TRUE),
      mu_bar_per_callable_sem = sem(r_callab),
      # New convention (per reference bp): biologically meaningful per-bp rate
      mu_bar_per_total        = mean(r_total, na.rm = TRUE),
      mu_bar_per_total_sem    = sem(r_total)
    )
  }))
}))

print(simple_out, digits = 4)
write.csv(simple_out, "rate_table_simple.csv", row.names = FALSE)
cat("Saved: rate_table_simple.csv\n")

# =============================================================================
# STEP 4: BOOTSTRAP WITH GENOME-WIDE FtR CORRECTION (mu_hat)
# Generative model:    m_obs ~ Poisson(mu × L_total × t × (1 - p_ftr))
# Estimator:           mu_hat = m_corrected / (L_total × t)
#                              where m_corrected = m + Pois(m × p/(1-p))
# For each bootstrap replicate:
#   1. Resample lines with replacement (within strain)
#   2. Apply FtR correction to counts (numerator only; uses corrected lambda)
#   3. Compute rate using L_total × generations (NOT L_callable)
# =============================================================================
cat("\n--- Running bootstrap (N =", N_BOOT, ", genome-wide FtR) ---\n")
cat("    Denominator: L_total × generations (NOT L_callable × generations)\n")
cat("    Correction:  m_hat = m + Pois(m × p_ftr / (1 - p_ftr)), E ~ m/(1-p)\n\n")

run_bootstrap <- function(df, strain_name) {
  n       <- nrow(df)
  p_snp   <- p_ftr[[paste0(strain_name, "_SNP")]]
  p_indel <- p_ftr[[paste0(strain_name, "_Indel")]]

  boot_mat <- matrix(NA_real_, nrow = N_BOOT, ncol = 4,
                     dimnames = list(NULL, mut_classes))

  for (b in seq_len(N_BOOT)) {
    idx  <- sample.int(n, size = n, replace = TRUE)
    samp <- df[idx, ]

    # Step 1: FtR correction on counts (numerator only, genome-wide formulation)
    c_snv   <- ftr_correct(samp$count_SNV,   p_snp)
    c_del   <- ftr_correct(samp$count_DEL,   p_indel)
    c_ins   <- ftr_correct(samp$count_INS,   p_indel)
    c_total <- c_snv + c_del + c_ins

    # Step 2: rate = corrected_count / (L_total × generations)
    # L_total is a CONSTANT (the reference genome length), NOT line-specific callable
    den <- denom_total(samp$gens)
    boot_mat[b, "SNV"]   <- mean(c_snv   / den, na.rm = TRUE)
    boot_mat[b, "DEL"]   <- mean(c_del   / den, na.rm = TRUE)
    boot_mat[b, "INS"]   <- mean(c_ins   / den, na.rm = TRUE)
    boot_mat[b, "TOTAL"] <- mean(c_total / den, na.rm = TRUE)
  }
  boot_mat
}

boot_HK <- run_bootstrap(mt %>% filter(strain == "HK104"), "HK104")
boot_PB <- run_bootstrap(mt %>% filter(strain == "PB800"), "PB800")
cat("Bootstrap complete.\n")

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
# STEP 5: COMBINED DISPLAY TABLE
# Rates reported × 10^9
# =============================================================================
cat("\n--- Building display table ---\n")
scale <- 1e9

display <- bind_rows(lapply(c("HK104", "PB800"), function(str) {
  s_row <- simple_out  %>% filter(strain == str)
  b_row <- boot_summary %>% filter(strain == str)
  meta  <- s_row %>% filter(class == "SNV") %>%
    select(N_lines, pct_genome, pct_genome_sem, mean_gens, gens_sem)

  bind_rows(lapply(mut_classes, function(cl) {
    s <- s_row  %>% filter(class == cl)
    b <- b_row  %>% filter(class == cl)
    data.frame(
      Strain             = str,
      Class              = cl,
      N_Lines            = meta$N_lines,
      pct_Genome         = round(meta$pct_genome,     3),
      pct_Genome_SEM     = round(meta$pct_genome_sem, 3),
      t_bar              = round(meta$mean_gens,      1),
      t_bar_SEM          = round(meta$gens_sem,       2),
      # Old convention (per-callable, NO correction)
      mu_bar_callab_x1e9 = round(s$mu_bar_per_callable     * scale, 4),
      mu_bar_callab_SEM  = round(s$mu_bar_per_callable_sem * scale, 4),
      # New convention (per-reference-bp, NO correction)
      mu_bar_total_x1e9  = round(s$mu_bar_per_total        * scale, 4),
      mu_bar_total_SEM   = round(s$mu_bar_per_total_sem    * scale, 4),
      # FtR-corrected mu_hat (per-reference-bp, genome-wide FtR model)
      mu_hat_x1e9        = round(b$mu_hat * scale, 4),
      mu_hat_CI_lo       = round(b$ci_lo  * scale, 4),
      mu_hat_CI_hi       = round(b$ci_hi  * scale, 4)
    )
  }))
}))

print(display, digits = 4)
write.csv(display, "rate_table_display.csv", row.names = FALSE)
cat("Saved: rate_table_display.csv\n")

# =============================================================================
# SUMMARY -- model and expected numerical shift vs. v1
# =============================================================================
cat("\n========================================================================\n")
cat("SUMMARY -- expected shift from v1 (within-callable FtR) to v2 (genome-wide)\n")
cat("========================================================================\n")
cat("v1 (old): mu_hat_old = m × (1 + p_ftr) / (L_callable × t)\n")
cat("v2 (new): mu_hat_new = m / [(1 - p_ftr) × L_total × t]\n\n")
cat("With L_callable/L_total ≈ 0.91 and p_ftr ≈ 0.05:\n")
cat("    mu_hat_new / mu_hat_old  ≈  (L_callable/L_total) × (1+p) × (1-p)\n")
cat("                              =  (L_callable/L_total) × (1 - p^2)\n")
cat("                              ≈  0.91 × 0.997\n")
cat("                              ≈  0.91\n\n")
cat("So v2 mu_hat values are ~9% LOWER than v1 in absolute numerical value,\n")
cat("but represent the per-bp mutation rate of the REFERENCE GENOME -- the\n")
cat("biologically meaningful quantity for cross-species comparison.\n\n")
cat("The HK104:PB800 ratios are essentially unchanged (within rounding) because\n")
cat("the correction factor is similar for both strains.\n\n")
cat("Cross-check: if v1 was reported as 'per callable site' and v2 as 'per\n")
cat("reference bp', they are different units of the same underlying mu.\n")
