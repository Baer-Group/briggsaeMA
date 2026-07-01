# =============================================================================
# §3.4.1 — HK104 vs PB800 per-line rate permutation tests
# C. briggsae MA experiment, 81-line cohort (32 HK104 + 49 PB800)
#
# v2 fix (2026-05-24): cast SR_callable_3x to numeric before multiplying with
# Generations. Without this cast, R's 32-bit integer arithmetic overflows
# (97 million × 238 ≈ 23 billion > .Machine$integer.max ≈ 2.1 billion) and
# silently returns NA, which then propagates into Del_rate, Ins_rate, and
# the subsequent permutation test outputs for those two classes.
#
# Test: difference of mean per-line rates between strains, under label
#       reshuffle. 10,000 permutations, two-sided p.
# Seed: 20260522 (same value as the Python run). Note: R's sample() and
#       NumPy's default_rng() use different RNG implementations, so p-values
#       will match to ~2 decimals but not bit-exactly. Point estimates
#       (means, ratios, delta) are deterministic and will match exactly.
# =============================================================================
suppressWarnings(try({
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}, silent = TRUE))
cat("Working directory:", getwd(), "\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

set.seed(20260522)
N_PERM <- 10000L
INPUT  <- "master_table.csv"          # adjust path if needed
OUTPUT <- "permutation_results_R.csv"

# ---- Load and prepare per-line rates -----------------------------------------
mt <- read.csv(INPUT, stringsAsFactors = FALSE)

# Cast columns used in rate denominators to numeric, preventing 32-bit
# integer overflow when callable_bp × generations exceeds ~2.1e9.
mt$SR_callable_3x  <- as.numeric(mt$SR_callable_3x)
mt$SR_callable_10x <- as.numeric(mt$SR_callable_10x)   # not used below but cast for safety
mt$Generations     <- as.numeric(mt$Generations)

stopifnot(nrow(mt) == 81L, all(mt$Generations == 238))

# DEL and INS aren't pre-computed as rates; derive from raw counts at DP >= 3.
del_cols <- c("SR_3x_del_1","SR_3x_del_2","SR_3x_del_3_5",
              "SR_3x_del_6_10","SR_3x_del_10plus")
ins_cols <- c("SR_3x_ins_1","SR_3x_ins_2","SR_3x_ins_3_5",
              "SR_3x_ins_6_10","SR_3x_ins_10plus")

# Defensive: ensure each count column is numeric too (CSV "NA" strings would
# otherwise make a column character, which silently breaks rowSums).
for (c in c(del_cols, ins_cols)) {
  mt[[c]] <- suppressWarnings(as.numeric(mt[[c]]))
}

mt$n_del    <- rowSums(mt[, del_cols], na.rm = TRUE)
mt$n_ins    <- rowSums(mt[, ins_cols], na.rm = TRUE)

# Per-line rates -- denominator is now numeric, so no overflow risk.
mt$Del_rate <- mt$n_del / (mt$SR_callable_3x * mt$Generations)
mt$Ins_rate <- mt$n_ins / (mt$SR_callable_3x * mt$Generations)

# Sanity check: cohort sizes and no NA in computed rate columns
cat(sprintf("Cohort: HK104 N=%d, PB800 N=%d, total N=%d\n",
            sum(mt$Strain == "HK104"),
            sum(mt$Strain == "PB800"),
            nrow(mt)))
stopifnot(!any(is.na(mt$Del_rate)),
          !any(is.na(mt$Ins_rate)),
          !any(is.na(mt$SR_SNP_rate_3x)),
          !any(is.na(mt$SR_Indel_rate_3x)),
          !any(is.na(mt$SR_Total_rate_3x)))
cat("Sanity check passed: no NA in any rate column.\n")

# ---- Permutation engine ------------------------------------------------------
# Test stat: mean(HK104 rates) - mean(PB800 rates).
# Two-sided p = fraction of |null_T| >= |T_obs|.
perm_test <- function(rates, strain, n_perm = N_PERM) {
  hk <- strain == "HK104"
  T_obs <- mean(rates[hk]) - mean(rates[!hk])
  null_T <- replicate(n_perm, {
    s <- sample(strain)              # reshuffle labels among 81 lines
    mean(rates[s == "HK104"]) - mean(rates[s == "PB800"])
  })
  p_two <- mean(abs(null_T) >= abs(T_obs))
  list(mean_HK = mean(rates[hk]),
       mean_PB = mean(rates[!hk]),
       T_obs   = T_obs,
       p_two   = p_two)
}

# ---- Run for all 5 classes ---------------------------------------------------
classes <- list(
  SNV       = "SR_SNP_rate_3x",
  Deletion  = "Del_rate",
  Insertion = "Ins_rate",
  Indel     = "SR_Indel_rate_3x",
  Total     = "SR_Total_rate_3x"
)

results <- lapply(names(classes), function(cls) {
  r <- perm_test(mt[[ classes[[cls]] ]], mt$Strain)
  data.frame(Class            = cls,
             mean_HK_x1e9     = r$mean_HK * 1e9,
             mean_PB_x1e9     = r$mean_PB * 1e9,
             HK_over_PB       = r$mean_HK / r$mean_PB,
             Delta_x1e9       = r$T_obs   * 1e9,
             p_two_sided      = r$p_two,
             p_reported       = ifelse(r$p_two == 0,
                                       sprintf("< %.4f", 1/N_PERM),
                                       sprintf("= %.4f", r$p_two)))
})
results <- do.call(rbind, results)

# ---- Print and save ----------------------------------------------------------
print(results, row.names = FALSE, digits = 4)
write.csv(results, OUTPUT, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", OUTPUT))

# Expected output with seed 20260522 and N_PERM = 10000 (point estimates
# match Python replication exactly; p-values will be close but not bit-exact
# because R's sample() and NumPy's default_rng() differ):
# -----------------------------------------------------------------------------
#      Class mean_HK_x1e9 mean_PB_x1e9 HK_over_PB Delta_x1e9 p_two_sided p_reported
#        SNV       1.6550       1.7735     0.9332    -0.1185      ~0.18    = ~0.18
#   Deletion       0.3209       0.2710     1.1840     0.0499      ~0.08    = ~0.08
#  Insertion       0.5021       0.2848     1.7628     0.2173      0.0000   < 0.0001
#      Indel       0.8228       0.5558     1.4805     0.2669      0.0000   < 0.0001
#      Total       2.4791       2.3292     1.0644     0.1499      ~0.20    = ~0.20
