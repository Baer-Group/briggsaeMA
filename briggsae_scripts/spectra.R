# =============================================================================
# Mutation Rate Spectra - C. briggsae MA Analysis
# =============================================================================
# Input:  master_table.csv (in same directory as this script)
# Output: 4 PDF figures (+ PNG versions)
#
#   Fig 1: DP >= 3x  -- left panel = SNP spectrum, right panel = Indel spectrum
#                       (HK104 vs PB800, SR all lines)
#   Fig 2: DP >= 10x -- same layout (HK104 vs PB800, SR all lines)
#   Fig 3: DP >= 3x  -- left = SNP, right = Indel (SR vs LR, matched lines)
#   Fig 4: DP >= 10x -- same layout (SR vs LR, matched lines)
#
# Rate normalisation:
#   SNP GC-context classes : count / (callable_sites x GC_FRAC x G)
#   SNP AT-context classes : count / (callable_sites x AT_FRAC x G)
#   Indel classes          : count / (callable_sites x G)
#
# Indel order (signed: negative = deletion, positive = insertion):
#   <-10 | -10 to -6 | -5 to -3 | -2 | -1 || 1 | 2 | 3 to 5 | 6 to 10 | 10<
#
# Error bars: SEM (standard error of the mean)
# =============================================================================

library(rstudioapi)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

setwd(dirname(getActiveDocumentContext()$path))
getwd()

# =============================================================================
# CONSTANTS
# =============================================================================
GENERATIONS <- 238
GC_FRAC     <- 0.37
AT_FRAC     <- 1 - GC_FRAC   # 0.63

# SNP classes and their genomic context for normalisation
SNP_CLASSES <- c("GC_to_AT", "GC_to_CG", "GC_to_TA", "AT_to_GC", "AT_to_CG", "AT_to_TA")
SNP_CONTEXT <- c("GC", "GC", "GC", "AT", "AT", "AT")
SNP_LABELS  <- c("GC->AT", "GC->CG", "GC->TA", "AT->GC", "AT->CG", "AT->TA")

# Indel bins -- signed format: negative = deletion, positive = insertion
INDEL_BINS   <- c("del_10plus", "del_6_10", "del_3_5", "del_2", "del_1",
                  "ins_1", "ins_2", "ins_3_5", "ins_6_10", "ins_10plus")
INDEL_LABELS <- c("<-10", "-10 to -6", "-5 to -3", "-2", "-1",
                  "1", "2", "3 to 5", "6 to 10", "10<")

# Colours
COL_HK     <- "blue4"
COL_PB     <- "coral"
COL_HK_LR <- "steelblue"
COL_PB_LR <- "lightsalmon"

# =============================================================================
# LOAD DATA
# =============================================================================
cat("Loading master table...\n")
mt <- read.csv("master_table.csv", stringsAsFactors = FALSE)
cat(sprintf("  %d lines loaded\n", nrow(mt)))

# =============================================================================
# HELPERS: compute per-line per-class rates -- long format
# =============================================================================
compute_snp_rates <- function(df, prefix, callable_col, platform) {
  out  <- df %>% filter(!is.na(.data[[callable_col]]), !is.na(Strain))
  rows <- list()
  for (i in seq_along(SNP_CLASSES)) {
    cl      <- SNP_CLASSES[i]
    frac    <- ifelse(SNP_CONTEXT[i] == "GC", GC_FRAC, AT_FRAC)
    col_cnt <- paste0(prefix, "_", cl)
    rows[[i]] <- out %>%
      select(Line, Strain,
             count    = all_of(col_cnt),
             callable = all_of(callable_col)) %>%
      mutate(rate     = count / (callable * frac * GENERATIONS),
             class    = cl,
             platform = platform)
  }
  bind_rows(rows) %>% select(Line, Strain, class, rate, platform)
}

compute_indel_rates <- function(df, prefix, callable_col, platform) {
  out  <- df %>% filter(!is.na(.data[[callable_col]]), !is.na(Strain))
  rows <- list()
  for (i in seq_along(INDEL_BINS)) {
    bn      <- INDEL_BINS[i]
    col_cnt <- paste0(prefix, "_", bn)
    rows[[i]] <- out %>%
      select(Line, Strain,
             count    = all_of(col_cnt),
             callable = all_of(callable_col)) %>%
      mutate(rate     = count / (callable * GENERATIONS),
             class    = bn,
             platform = platform)
  }
  bind_rows(rows) %>% select(Line, Strain, class, rate, platform)
}

# =============================================================================
# BUILD LONG DATA FRAMES
# =============================================================================
cat("Computing per-line rates...\n")

# SR all lines
sr_snp_3x  <- compute_snp_rates(mt,   "SR_3x",  "SR_callable_3x",  "SR")
sr_snp_10x <- compute_snp_rates(mt,   "SR_10x", "SR_callable_10x", "SR")
sr_ind_3x  <- compute_indel_rates(mt, "SR_3x",  "SR_callable_3x",  "SR")
sr_ind_10x <- compute_indel_rates(mt, "SR_10x", "SR_callable_10x", "SR")

# Matched lines only (have both SR and LR)
mt_m <- mt %>% filter(!is.na(LR_sample))

lr_snp_3x  <- bind_rows(
  compute_snp_rates(mt_m, "SR_matched_3x",  "SR_callable_3x",  "SR"),
  compute_snp_rates(mt_m, "LR_3x",          "LR_callable_3x",  "LR"))
lr_snp_10x <- bind_rows(
  compute_snp_rates(mt_m, "SR_matched_10x", "SR_callable_10x", "SR"),
  compute_snp_rates(mt_m, "LR_10x",         "LR_callable_10x", "LR"))
lr_ind_3x  <- bind_rows(
  compute_indel_rates(mt_m, "SR_matched_3x",  "SR_callable_3x",  "SR"),
  compute_indel_rates(mt_m, "LR_3x",          "LR_callable_3x",  "LR"))
lr_ind_10x <- bind_rows(
  compute_indel_rates(mt_m, "SR_matched_10x", "SR_callable_10x", "SR"),
  compute_indel_rates(mt_m, "LR_10x",         "LR_callable_10x", "LR"))

# =============================================================================
# SUMMARISE: mean +/- SEM per group
# =============================================================================
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

summarise_strain <- function(df) {
  df %>%
    group_by(Strain, class) %>%
    summarise(mean_rate = mean(rate, na.rm = TRUE),
              se_rate   = sem(rate),
              n         = sum(!is.na(rate)),
              .groups   = "drop")
}

summarise_platform <- function(df) {
  df %>%
    mutate(group = factor(paste0(Strain, "-", platform),
                          levels = c("HK104-SR", "HK104-LR", "PB800-SR", "PB800-LR"))) %>%
    group_by(group, class) %>%
    summarise(mean_rate = mean(rate, na.rm = TRUE),
              se_rate   = sem(rate),
              n         = sum(!is.na(rate)),
              .groups   = "drop")
}

# Apply factor levels
set_snp_levels <- function(df) {
  df %>% mutate(class = factor(class, levels = SNP_CLASSES, labels = SNP_LABELS))
}
set_ind_levels <- function(df) {
  df %>% mutate(class = factor(class, levels = INDEL_BINS, labels = INDEL_LABELS))
}

ss3_snp  <- summarise_strain(sr_snp_3x)  %>% set_snp_levels()
ss10_snp <- summarise_strain(sr_snp_10x) %>% set_snp_levels()
ss3_ind  <- summarise_strain(sr_ind_3x)  %>% set_ind_levels()
ss10_ind <- summarise_strain(sr_ind_10x) %>% set_ind_levels()

lp3_snp  <- summarise_platform(lr_snp_3x)  %>% set_snp_levels()
lp10_snp <- summarise_platform(lr_snp_10x) %>% set_snp_levels()
lp3_ind  <- summarise_platform(lr_ind_3x)  %>% set_ind_levels()
lp10_ind <- summarise_platform(lr_ind_10x) %>% set_ind_levels()

# =============================================================================
# SHARED THEME
# =============================================================================
theme_spectra <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y  = element_text(size = 9),
      axis.title   = element_text(size = 10),
      legend.position = "top",
      legend.title    = element_text(size = 9),
      legend.text     = element_text(size = 9),
      plot.title      = element_blank()
    )
}

# =============================================================================
# PANEL BUILDERS
# =============================================================================

# SNP panel -- strain comparison (HK104 vs PB800)
snp_strain_panel <- function(df) {
  ggplot(df, aes(x = class, y = mean_rate, fill = Strain)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.7),
             width = 0.65, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_rate - se_rate,
                      ymax = mean_rate + se_rate),
                  position = position_dodge(width = 0.7),
                  width = 0.22, linewidth = 0.5) +
    scale_fill_manual(values = c("HK104" = COL_HK, "PB800" = COL_PB),
                      name = "Strain") +
    scale_y_continuous(labels = scales::scientific,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Substitution class",
         y = "Rate (per site per gen)") +
    theme_spectra()
}

# Indel panel -- strain comparison (HK104 vs PB800)
ind_strain_panel <- function(df) {
  ggplot(df, aes(x = class, y = mean_rate, fill = Strain)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.7),
             width = 0.65, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_rate - se_rate,
                      ymax = mean_rate + se_rate),
                  position = position_dodge(width = 0.7),
                  width = 0.22, linewidth = 0.5) +
    geom_vline(xintercept = 5.5, linetype = "dashed",
               color = "grey50", linewidth = 0.6) +
    scale_fill_manual(values = c("HK104" = COL_HK, "PB800" = COL_PB),
                      name = "Strain") +
    scale_y_continuous(labels = scales::scientific,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Indel class",
         y = "Rate (per site per gen)") +
    theme_spectra()
}

# SNP panel -- platform comparison (SR vs LR, matched lines)
snp_platform_panel <- function(df) {
  ggplot(df, aes(x = class, y = mean_rate, fill = group)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8),
             width = 0.75, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_rate - se_rate,
                      ymax = mean_rate + se_rate),
                  position = position_dodge(width = 0.8),
                  width = 0.22, linewidth = 0.5) +
    scale_fill_manual(
      values = c("HK104-SR" = COL_HK,    "HK104-LR" = COL_HK_LR,
                 "PB800-SR" = COL_PB,    "PB800-LR" = COL_PB_LR),
      name = "Strain - Platform") +
    scale_y_continuous(labels = scales::scientific,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Substitution class",
         y = "Rate (per site per gen)") +
    theme_spectra()
}

# Indel panel -- platform comparison (SR vs LR, matched lines)
ind_platform_panel <- function(df) {
  ggplot(df, aes(x = class, y = mean_rate, fill = group)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8),
             width = 0.75, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_rate - se_rate,
                      ymax = mean_rate + se_rate),
                  position = position_dodge(width = 0.8),
                  width = 0.22, linewidth = 0.5) +
    geom_vline(xintercept = 5.5, linetype = "dashed",
               color = "grey50", linewidth = 0.6) +
    scale_fill_manual(
      values = c("HK104-SR" = COL_HK,    "HK104-LR" = COL_HK_LR,
                 "PB800-SR" = COL_PB,    "PB800-LR" = COL_PB_LR),
      name = "Strain - Platform") +
    scale_y_continuous(labels = scales::scientific,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Indel class",
         y = "Rate (per site per gen)") +
    theme_spectra()
}

# =============================================================================
# ASSEMBLE FIGURES (SNP | Indel side by side via patchwork)
# =============================================================================
cat("Assembling figures...\n")

# Fig 1: DP >= 3x, HK104 vs PB800 (SR all lines)
fig1 <- snp_strain_panel(ss3_snp) +
        ind_strain_panel(ss3_ind) +
        plot_layout(guides = "collect") &
        theme(legend.position = "top")

# Fig 2: DP >= 10x, HK104 vs PB800 (SR all lines)
fig2 <- snp_strain_panel(ss10_snp) +
        ind_strain_panel(ss10_ind) +
        plot_layout(guides = "collect") &
        theme(legend.position = "top")

# Fig 3: DP >= 3x, SR vs LR (matched lines)
fig3 <- snp_platform_panel(lp3_snp) +
        ind_platform_panel(lp3_ind) +
        plot_layout(guides = "collect") &
        theme(legend.position = "top")

# Fig 4: DP >= 10x, SR vs LR (matched lines)
fig4 <- snp_platform_panel(lp10_snp) +
        ind_platform_panel(lp10_ind) +
        plot_layout(guides = "collect") &
        theme(legend.position = "top")

# =============================================================================
# SAVE FIGURES
# =============================================================================
cat("Saving figures...\n")

save_fig <- function(p, name, w = 12, h = 5) {
  ggsave(paste0(name, ".pdf"), plot = p, width = w, height = h)
  ggsave(paste0(name, ".png"), plot = p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s.pdf / .png\n", name))
}

save_fig(fig1, "Fig1_Spectra_SR_DP3x_HKvsPB")
save_fig(fig2, "Fig2_Spectra_SR_DP10x_HKvsPB")
save_fig(fig3, "Fig3_Spectra_SR_vs_LR_DP3x",  w = 13, h = 5)
save_fig(fig4, "Fig4_Spectra_SR_vs_LR_DP10x", w = 13, h = 5)

cat("\nDone. 4 figures saved (PDF + PNG).\n")
