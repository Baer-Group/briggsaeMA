# =============================================================================
# Mononucleotide Repeat Analysis -- C. briggsae MA Lines
# Framework: Rajaei et al. 2021 (Genome Research)
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
#   Fig_SNV_rates.pdf/.png    -- Fig 2A-C equiv: SNV rates by context
#   Fig_Indel_rates.pdf/.png  -- Fig 2D-F equiv: indel rates by context
#   Fig_SNV_spectrum.pdf/.png -- Fig 4A-C equiv: proportions by context
#   Fig_rate_ratio.pdf/.png   -- summary: rate ratio with bootstrap CI
#   Fig_TsTv_context.pdf/.png
#   Fig_unit_heatmap.pdf/.png -- EXTENSION: per-unit enrichment
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
REPEAT_BP   <- 6001892
P_MONO      <- REPEAT_BP / GENOME_SIZE   # 0.05652: fraction of genome in mono
P_NONMONO   <- 1 - P_MONO               # 0.94348

GC_FRAC     <- 0.37
AT_FRAC     <- 1 - GC_FRAC              # 0.63
GENERATIONS <- 238
N_BOOT      <- 10000
set.seed(2025)

SHEET       <- "3x"

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

theme_rajaei <- function(base_size=11) {
  theme_bw(base_size=base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x  = element_text(angle=45, hjust=1, size=9),
      axis.text.y  = element_text(size=9),
      axis.title   = element_text(size=10),
      legend.position  = "top",
      legend.title     = element_text(size=9),
      legend.text      = element_text(size=9),
      strip.background = element_rect(fill="grey92"),
      strip.text       = element_text(size=9, face="bold"),
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

bed <- fread("briggsae_mono_repeats_with_unit.bed",
             col.names=c("CHROM","start0","end1","repeat_unit"),
             header=FALSE)
bed[, repeat_unit := toupper(trimws(repeat_unit))]
bed[, start1 := start0 + 1]
setkey(bed, CHROM, start1, end1)

# Per-unit bp for extension analysis
unit_bp <- bed[, .(unit_bp = sum(end1-start0)), by=repeat_unit]
unit_bp[, p_null_unit := unit_bp / GENOME_SIZE]
cat("Per-unit null:\n"); print(unit_bp)

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

load_mutations <- function(file) {
  df <- as.data.table(read_excel(file, sheet=SHEET))

  # Strain
  df[, sample_num := as.integer(sub("_.*","", as.character(Sample)))]
  df[, Strain := ifelse(sample_num <= 298, "HK104", "PB800")]
  df[, Line   := as.character(sample_num)]
  df[, sample_num := NULL]

  # foverlaps for repeat_unit
  df[, pos_end := POS]
  setkey(df, CHROM, POS, pos_end)
  hits <- foverlaps(df,
                    bed[, .(CHROM, start1, end1, repeat_unit)],
                    by.x=c("CHROM","POS","pos_end"),
                    by.y=c("CHROM","start1","end1"),
                    type="within", nomatch=NA)
  df[, repeat_unit := hits$repeat_unit]
  df[, pos_end := NULL]

  # Classify mutation class
  df <- df %>%
    rowwise() %>%
    mutate(mut_class = case_when(
      TYPE=="snp"   ~ classify_snv(REF, ALT),
      TYPE=="indel" ~ classify_indel(REF, ALT),
      TRUE          ~ NA_character_
    )) %>%
    ungroup()

  # Ensure mono_repeat is logical
  df$mono_repeat <- as.logical(df$mono_repeat)
  df
}

sr    <- load_mutations("SR_final_mutation_list_mono_classified.xlsx")
sr_lr <- load_mutations("SR_LR_final_mutation_list_mono_classified.xlsx")

cat(sprintf("  SR:    %d mutations | %d in mono repeats\n",
            nrow(sr), sum(sr$mono_repeat, na.rm=TRUE)))
cat(sprintf("  SR+LR: %d mutations | %d in mono repeats\n",
            nrow(sr_lr), sum(sr_lr$mono_repeat, na.rm=TRUE)))

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
      # GC/AT fraction for normalisation
      frac = case_when(
        mut_class %in% names(SNP_FRAC) ~ SNP_FRAC[mut_class],
        TRUE                           ~ 1.0   # indels: no GC/AT normalisation
      ),
      # Context-specific callable sites
      callable_ctx = case_when(
        context == "mono"    ~ callable * P_MONO,
        context == "nonmono" ~ callable * P_NONMONO,
        context == "total"   ~ as.numeric(callable)
      ),
      # Rate per site per generation x10^9
      rate = count / (callable_ctx * frac * GENERATIONS)
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
write.csv(rate_table, "rate_table_rajaei.csv", row.names=FALSE)
cat("Saved: rate_table_rajaei.csv\n")

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

write.csv(ratio_results, "rate_ratio_bootstrap.csv", row.names=FALSE)
cat("Saved: rate_ratio_bootstrap.csv\n")

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
write.csv(spectrum_tests, "spectrum_stats.csv", row.names=FALSE)
cat("Saved: spectrum_stats.csv\n")

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
  plot_layout(guides="collect") & theme(legend.position="top")

save_fig(fig_snv, "Fig_SNV_rates", 14, 5)

# -- Fig 2: Indel rates by context (Rajaei Fig 2D-F equivalent) -------------
fig_ind <- make_rate_figure(INDEL_BINS, INDEL_LABELS, add_vline=TRUE) +
  plot_layout(guides="collect") & theme(legend.position="top")

save_fig(fig_ind, "Fig_Indel_rates", 16, 5)

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
  plot_layout(guides="collect") & theme(legend.position="top")

save_fig(fig_spec, "Fig_SNV_spectrum", 14, 5)

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
    theme(axis.text.x=element_text(angle=45, hjust=1))

save_fig(fig_ratio, "Fig_rate_ratio", 13, 5)

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

save_fig(fig_tstv, "Fig_TsTv_context", 6, 4)

# =============================================================================
# EXTENSION: Per-unit enrichment heatmap (beyond Rajaei)
# For each SNV class x repeat unit (A/T/C/G):
#   null    = unit_bp / genome_bp  (unit-specific)
#   obs     = N(class in unit) / N(class total)
#   ratio   = obs / null
#   test    = binomial, BH-corrected
# =============================================================================
cat("\n--- Extension: Per-unit enrichment (beyond Rajaei) ---\n")

unit_null_vec <- setNames(unit_bp$p_null_unit, unit_bp$repeat_unit)

unit_enr <- bind_rows(lapply(c("HK104","PB800"), function(str) {
  snv <- sr %>%
    filter(Strain==str, TYPE=="snp",
           !is.na(mut_class), mut_class %in% SNP_CLASSES)

  bind_rows(lapply(SNP_CLASSES, function(cl) {
    sub   <- snv %>% filter(mut_class==cl)
    n_tot <- nrow(sub)
    bind_rows(lapply(UNIT_LEVELS, function(unit) {
      p_u  <- unit_null_vec[unit]
      if (is.na(p_u) || n_tot==0) return(NULL)
      n_in <- sum(sub$repeat_unit==unit, na.rm=TRUE)
      bt   <- binom.test(n_in, n_tot, p=p_u, alternative="two.sided")
      data.frame(
        Strain=str, class=cl,
        label=SNP_LABELS[match(cl,SNP_CLASSES)],
        unit=unit, n_total=n_tot, n_in_unit=n_in,
        obs_frac=n_in/n_tot, null_frac=p_u,
        enrichment_ratio=(n_in/n_tot)/p_u,
        p_value=bt$p.value, p_adj=NA, significant=FALSE
      )
    }))
  }))
})) %>%
  group_by(Strain) %>%
  mutate(p_adj=p.adjust(p_value, method="BH"),
         significant=!is.na(p_adj) & p_adj < 0.05) %>%
  ungroup()

write.csv(unit_enr, "unit_enrichment_extension.csv", row.names=FALSE)

fig_heatmap <- unit_enr %>%
  mutate(
    label    = factor(label, levels=SNP_LABELS),
    unit     = factor(unit,  levels=UNIT_LEVELS),
    log2_enr = log2(enrichment_ratio),
    sig_lbl  = ifelse(significant, "*", "")
  ) %>%
  ggplot(aes(x=unit, y=label, fill=log2_enr)) +
    geom_tile(color="white", linewidth=0.4) +
    geom_text(aes(label=sig_lbl), size=5, vjust=0.8) +
    scale_fill_gradient2(
      low="#1F77B4", mid="white", high="#D62728", midpoint=0,
      name="log2(enrichment)\n* FDR<0.05"
    ) +
    facet_wrap(~Strain, nrow=1) +
    labs(x="Repeat unit", y="SNV class") +
    theme_rajaei() +
    theme(axis.text.x=element_text(angle=0, hjust=0.5),
          panel.grid=element_blank())

save_fig(fig_heatmap, "Fig_unit_heatmap", 9, 5)

# =============================================================================
# DONE
# =============================================================================
cat("\n=== Complete. All outputs saved. ===\n")
cat("\n-- Rajaei-framework outputs: --\n")
cat("  rate_table_rajaei.csv      (Table 1 equivalent)\n")
cat("  rate_ratio_bootstrap.csv   (mu_mono/mu_nonmono + 95% CI)\n")
cat("  spectrum_stats.csv         (MC-FET p-values)\n")
cat("  Fig_SNV_rates.pdf/.png     (Fig 2A-C equivalent)\n")
cat("  Fig_Indel_rates.pdf/.png   (Fig 2D-F equivalent)\n")
cat("  Fig_SNV_spectrum.pdf/.png  (Fig 4A-C equivalent)\n")
cat("  Fig_rate_ratio.pdf/.png    (summary: ratio + CI)\n")
cat("  Fig_TsTv_context.pdf/.png\n")
cat("\n-- Extension output: --\n")
cat("  unit_enrichment_extension.csv\n")
cat("  Fig_unit_heatmap.pdf/.png  (A/T/C/G run context)\n")

cat("\n-- Key numbers to compare with Rajaei Table 1 (C. elegans N2): --\n")
cat("  AT->TA nonmono: ~0.44 x10^-9 | mono: ~3.34 x10^-9 | ratio: ~7.6x\n")
cat("  Del-1  nonmono: ~0.19 x10^-9 | mono: ~1.43 x10^-9 | ratio: ~7.5x\n")
cat("  Ins-1  nonmono: ~0.07 x10^-9 | mono: ~1.07 x10^-9 | ratio: ~15x\n")
cat("  Your briggsae values are in rate_ratio_bootstrap.csv\n")
