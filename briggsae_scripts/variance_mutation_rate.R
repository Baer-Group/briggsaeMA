# =============================================================================
# Variance in Mutation Rate -- C. briggsae MA Lines
# HK104 vs PB800 (no O1/O2 hierarchy -- single-tier experiment)
#
# Framework: Saxena et al. 2019 / Saber et al. 2025 (sd_separate_o1o2.R)
# Adaptation: O1/O2 replaced by HK104/PB800 strain grouping
#
# Null hypothesis: all lines within a strain share one true mutation rate mu.
#   Simulate 10,000 datasets under this null (Poisson sampling + coverage
#   thinning), compute the SD of per-line rates in each replicate, and ask:
#   is the observed SD greater than the null expectation?
#   P-value = proportion of simulated SDs exceeding the observed SD.
#
# No FP correction: false positives handled upstream via replicate concordance
#   and IGV validation. FtR already corrected in rate_table_briggsae.R.
#   Simulation here uses raw observed counts (mu estimated from raw counts).
#
# Inputs:
#   master_table.csv
#
# Outputs:
#   variance_results.csv
#   Fig_variance_SNV.pdf/.png
#   Fig_variance_Indel.pdf/.png
#   Fig_variance_Total.pdf/.png
#   Fig_variance_combined.pdf/.png
# =============================================================================

library(rstudioapi)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

setwd(dirname(getActiveDocumentContext()$path))
set.seed(2025)

# =============================================================================
# CONSTANTS
# =============================================================================
GENOME_SIZE <- 106196309   # C. briggsae Feb2020 reference (bp)
GENERATIONS <- 238         # constant for all lines
NSIM        <- 10000

GC_FRAC <- 0.37
AT_FRAC <- 0.63

# master_table column names
COL_STRAIN   <- "Strain"
COL_CALLABLE <- "SR_callable_3x"

SNP_COLS  <- c("SR_3x_GC_to_AT","SR_3x_GC_to_CG","SR_3x_GC_to_TA",
               "SR_3x_AT_to_GC","SR_3x_AT_to_CG","SR_3x_AT_to_TA")
DEL_COLS  <- c("SR_3x_del_10plus","SR_3x_del_6_10","SR_3x_del_3_5",
               "SR_3x_del_2","SR_3x_del_1")
INS_COLS  <- c("SR_3x_ins_1","SR_3x_ins_2","SR_3x_ins_3_5",
               "SR_3x_ins_6_10","SR_3x_ins_10plus")

# Strain colors (matching dissertation palette)
COL_HK <- "blue4"
COL_PB <- "coral"

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================
cat("Loading master table...\n")
mt <- read.csv("master_table.csv", stringsAsFactors = FALSE)

# Derive strain if not present
if (!COL_STRAIN %in% names(mt)) {
  mt$sample_num <- as.integer(sub("_.*", "", as.character(mt$Line)))
  mt$Strain     <- ifelse(mt$sample_num <= 298, "HK104", "PB800")
}

# Aggregate counts per line
mt <- mt %>%
  filter(!is.na(.data[[COL_CALLABLE]]), !is.na(Strain)) %>%
  mutate(
    count_SNV   = rowSums(across(all_of(SNP_COLS)),  na.rm = TRUE),
    count_Indel = rowSums(across(all_of(c(DEL_COLS, INS_COLS))), na.rm = TRUE),
    count_Total = count_SNV + count_Indel,
    callable    = .data[[COL_CALLABLE]],
    cov_frac    = callable / GENOME_SIZE    # fraction of genome callable
  )

# Split by strain
hk <- mt %>% filter(Strain == "HK104")
pb <- mt %>% filter(Strain == "PB800")
cat(sprintf("  HK104: %d lines | PB800: %d lines\n", nrow(hk), nrow(pb)))

# =============================================================================
# HELPER: Observed per-line rate (per bp per gen) and SD
# =============================================================================
per_line_rates <- function(df, count_col) {
  df[[count_col]] / (df$callable * GENERATIONS)
}

observed_sd <- function(df, count_col) {
  sd(per_line_rates(df, count_col), na.rm = TRUE)
}

observed_mean <- function(df, count_col) {
  mean(per_line_rates(df, count_col), na.rm = TRUE)
}

# =============================================================================
# MLE FOR MUTATION RATE (ratio of sums -- pooled MLE)
# mu_genome_per_gen = total mutations / total generations
# mu_bp_per_gen     = total mutations / (total callable * total generations)
# Simulation uses mu_genome (per-genome per-generation) then thins by coverage
# =============================================================================
mle_genome <- function(df, count_col) {
  # MLE = total mutations / (n_lines x generations)
  # sum(GENERATIONS) would return the scalar 238, NOT n x 238 -- hence
  # we explicitly multiply nrow(df) x GENERATIONS as the denominator
  sum(df[[count_col]], na.rm = TRUE) / (nrow(df) * GENERATIONS)
}

# =============================================================================
# SIMULATION UNDER NULL (uniform rate within strain)
#
# For each replicate k:
#   For each line i:
#     1. Draw true mutations: raw ~ Poisson(GENERATIONS * mu_genome)
#        (mu_genome is per-genome per-generation)
#     2. Thin by coverage fraction: kept ~ Binomial(raw, cov_frac_i)
#        (only callable sites can harbour detected mutations)
#     3. Convert to rate: rate_i = kept / (cov_frac_i * GENERATIONS * GENOME_SIZE)
#   Compute SD of rates across lines
#
# NOTE: cov_frac and callable_sites serve different roles --
#   cov_frac enters thinning (step 2) and the denominator (step 3) separately.
#   They are NOT used simultaneously; step 2 is the numerator correction
#   and step 3 is the denominator normalisation.
# =============================================================================
simulate_null <- function(df, count_col, nsim = NSIM) {
  mu_g   <- mle_genome(df, count_col)   # per-genome per-generation
  n      <- nrow(df)
  cf     <- df$cov_frac                 # coverage fraction per line

  SD_sim   <- numeric(nsim)
  Mean_sim <- numeric(nsim)

  for (k in seq_len(nsim)) {
    rates <- numeric(n)
    for (i in seq_len(n)) {
      # Step 1: Poisson draw of true mutations per genome per experiment
      raw   <- rpois(1, GENERATIONS * mu_g)
      # Step 2: Thin by callable fraction (coverage)
      kept  <- rbinom(1, size = raw, prob = cf[i])
      # Step 3: Convert to per-bp per-gen rate
      rates[i] <- kept / (cf[i] * GENERATIONS * GENOME_SIZE)
    }
    SD_sim[k]   <- sd(rates,   na.rm = TRUE)
    Mean_sim[k] <- mean(rates, na.rm = TRUE)
  }
  list(sd = SD_sim, mean = Mean_sim)
}

# =============================================================================
# RUN SIMULATIONS (3 mutation classes x 2 strains = 6 runs)
# =============================================================================
cat("Running simulations (NSIM =", NSIM, ")...\n")
cat("  HK104 SNV...  ")
hk_snv  <- simulate_null(hk, "count_SNV");   cat("done\n")
cat("  HK104 Indel...")
hk_ind  <- simulate_null(hk, "count_Indel"); cat("done\n")
cat("  HK104 Total...")
hk_tot  <- simulate_null(hk, "count_Total"); cat("done\n")
cat("  PB800 SNV...  ")
pb_snv  <- simulate_null(pb, "count_SNV");   cat("done\n")
cat("  PB800 Indel...")
pb_ind  <- simulate_null(pb, "count_Indel"); cat("done\n")
cat("  PB800 Total...")
pb_tot  <- simulate_null(pb, "count_Total"); cat("done\n")

# =============================================================================
# OBSERVED SDs AND P-VALUES
# P = proportion of simulated SDs >= observed SD (right-tail test)
# =============================================================================
obs <- data.frame(
  class   = rep(c("SNV","Indel","Total"), 2),
  strain  = rep(c("HK104","PB800"), each = 3),
  obs_sd  = c(
    observed_sd(hk, "count_SNV"),   observed_sd(hk, "count_Indel"),
    observed_sd(hk, "count_Total"),
    observed_sd(pb, "count_SNV"),   observed_sd(pb, "count_Indel"),
    observed_sd(pb, "count_Total")
  ),
  obs_mean = c(
    observed_mean(hk, "count_SNV"),   observed_mean(hk, "count_Indel"),
    observed_mean(hk, "count_Total"),
    observed_mean(pb, "count_SNV"),   observed_mean(pb, "count_Indel"),
    observed_mean(pb, "count_Total")
  ),
  p_value = c(
    mean(hk_snv$sd >= observed_sd(hk, "count_SNV")),
    mean(hk_ind$sd >= observed_sd(hk, "count_Indel")),
    mean(hk_tot$sd >= observed_sd(hk, "count_Total")),
    mean(pb_snv$sd >= observed_sd(pb, "count_SNV")),
    mean(pb_ind$sd >= observed_sd(pb, "count_Indel")),
    mean(pb_tot$sd >= observed_sd(pb, "count_Total"))
  )
)

obs$mean_x1e9 <- obs$obs_mean * 1e9
obs$sd_x1e9   <- obs$obs_sd   * 1e9
obs$significant <- obs$p_value < 0.05

cat("\n=== Variance Results ===\n")
print(obs[, c("strain","class","obs_mean","obs_sd","p_value","significant")])
write.csv(obs, "variance_results.csv", row.names = FALSE)
cat("Saved: variance_results.csv\n")

# =============================================================================
# FIGURES: Overlay histogram of simulated SDs with observed values marked
# Structure mirrors sd_separate_o1o2.R: HK104 vs PB800 replaces O1 vs O2
# =============================================================================
plot_variance <- function(sim_hk, sim_pb, obs_hk, obs_pb,
                          p_hk, p_pb, title, file_out, bins = 60) {

  plot_data <- data.frame(
    sd     = c(sim_hk$sd, sim_pb$sd),
    strain = rep(c("HK104","PB800"), each = length(sim_hk$sd))
  )

  # significance label helper
  sig_label <- function(p) {
    if      (p < 0.001) "p < 0.001"
    else if (p < 0.01)  sprintf("p = %.3f", p)
    else if (p < 0.05)  sprintf("p = %.3f", p)
    else                sprintf("p = %.3f", p)
  }

  p <- ggplot(plot_data, aes(x = sd * 1e9, fill = strain)) +
    geom_histogram(position = "identity", alpha = 0.55, bins = bins,
                   color = "white", linewidth = 0.2) +

    # HK104 observed line
    geom_vline(xintercept = obs_hk * 1e9, color = COL_HK,
               linetype = "dashed", linewidth = 1.1) +
    annotate("text", x = obs_hk * 1e9, y = Inf,
             label = "Observed\nHK104", color = COL_HK,
             angle = 90, vjust = -0.4, hjust = 1.1, size = 4.5) +
    annotate("text", x = obs_hk * 1e9, y = 0,
             label = sig_label(p_hk), color = "grey20",
             vjust = 1.5, hjust = 0.5, size = 4) +

    # PB800 observed line
    geom_vline(xintercept = obs_pb * 1e9, color = COL_PB,
               linetype = "dashed", linewidth = 1.1) +
    annotate("text", x = obs_pb * 1e9, y = Inf,
             label = "Observed\nPB800", color = COL_PB,
             angle = 90, vjust = -0.4, hjust = 1.1, size = 4.5) +
    annotate("text", x = obs_pb * 1e9, y = 0,
             label = sig_label(p_pb), color = "grey20",
             vjust = 1.5, hjust = 0.5, size = 4) +

    scale_fill_manual(
      name   = "Strain (null simulation)",
      values = c(HK104 = COL_HK, PB800 = COL_PB)
    ) +
    labs(
      x     = expression(SD~of~per-line~rate~(x10^{-9}~"per bp per gen")),
      y     = "Frequency (simulated replicates)",
      title = title
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.title   = element_text(size = 14),
      axis.text    = element_text(size = 13),
      legend.text  = element_text(size = 13),
      legend.title = element_text(size = 13),
      plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.position = "top"
    )

  ggsave(file_out,          plot = p, width = 8, height = 5, dpi = 300)
  ggsave(sub(".pdf",".png", file_out), plot = p, width = 8, height = 5, dpi = 300)
  cat(sprintf("  Saved: %s\n", file_out))
  invisible(p)
}

cat("\nGenerating figures...\n")

p_snv <- plot_variance(
  hk_snv, pb_snv,
  observed_sd(hk, "count_SNV"),   observed_sd(pb, "count_SNV"),
  obs$p_value[obs$strain=="HK104" & obs$class=="SNV"],
  obs$p_value[obs$strain=="PB800" & obs$class=="SNV"],
  "SNV Mutation Rate Variance",
  "Fig_variance_SNV.pdf"
)

p_ind <- plot_variance(
  hk_ind, pb_ind,
  observed_sd(hk, "count_Indel"),  observed_sd(pb, "count_Indel"),
  obs$p_value[obs$strain=="HK104" & obs$class=="Indel"],
  obs$p_value[obs$strain=="PB800" & obs$class=="Indel"],
  "Indel Mutation Rate Variance",
  "Fig_variance_Indel.pdf"
)

p_tot <- plot_variance(
  hk_tot, pb_tot,
  observed_sd(hk, "count_Total"),  observed_sd(pb, "count_Total"),
  obs$p_value[obs$strain=="HK104" & obs$class=="Total"],
  obs$p_value[obs$strain=="PB800" & obs$class=="Total"],
  "Total Mutation Rate Variance",
  "Fig_variance_Total.pdf"
)

# Combined 3-panel figure
combined <- arrangeGrob(p_snv, p_ind, p_tot, ncol = 1)
ggsave("Fig_variance_combined.pdf", combined, width = 8, height = 15)
ggsave("Fig_variance_combined.png", combined, width = 8, height = 15, dpi = 300)
cat("  Saved: Fig_variance_combined.pdf/.png\n")

# =============================================================================
# SUMMARY PRINT
# =============================================================================
cat("\n=== SUMMARY ===\n")
cat(sprintf("Genome size:  %d bp\nGenerations:  %d\nNSIM:         %d\n\n",
            GENOME_SIZE, GENERATIONS, NSIM))

for (i in seq_len(nrow(obs))) {
  cat(sprintf("%-5s %-6s | mean = %.3e | SD = %.3e | p = %.4f %s\n",
              obs$strain[i], obs$class[i],
              obs$obs_mean[i], obs$obs_sd[i],
              obs$p_value[i],
              ifelse(obs$significant[i], "<-- SIGNIFICANT", "")))
}

cat("\n=== DRIFT-BARRIER INTERPRETATION ===\n")
cat("If p < 0.05: observed variance exceeds Poisson null\n")
cat("  -> among-line heterogeneity in mutation rate\n")
cat("  -> consistent with accumulation of (anti)mutator alleles\n")
cat("If p >= 0.05: no evidence of excess variance\n")
cat("  -> mutation rate appears homogeneous within strain\n")
cat("\nCross-species reference (Saber et al. 2025, C. elegans):\n")
cat("  Compare p-values and SD magnitudes between briggsae and elegans\n")
cat("  to assess whether variance input rate is conserved across species.\n")

cat("\nDone.\n")

