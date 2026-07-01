# =============================================================================
# Master Table Builder — C. briggsae MA Mutation Rate & Spectra Analysis
# =============================================================================
# Inputs (all in same directory as this script):
#   - C_briggsae_MetaData.xlsx
#   - SR_final_mutation_list.xlsx
#   - SR_LR_final_mutation_list.xlsx
#
# Output:
#   - master_table.csv  (one row per MA line)
#
# Notes:
#   - Reference genome length: 106,196,309 bp
#   - Generations: 238 (uniform; update when per-line data available)
#   - LR HK_270 = SR line 290; LR PB_342 = SR line 343
#   - HK_204 excluded from LR (coverage too low: 6.7x)
#   - Ancestors excluded from all calculations
#   - Depth pairing: 3x mutations <-> 3x callable sites; 10x <-> 10x
# =============================================================================

library(rstudioapi)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

setwd(dirname(getActiveDocumentContext()$path))
getwd()

# =============================================================================
# CONSTANTS
# =============================================================================
REF_LENGTH <- 106196309
GENERATIONS <- 238

# LR lines to exclude entirely
LR_EXCLUDE <- c("HK_204", "HK104_ANC_1", "PB800_ANC__1")

# LR sample ID -> SR line ID mapping
# (LR naming uses strain prefix; SR uses numeric only)
LR_TO_SR <- c(
  "HK_205" = "205", "HK_209" = "209", "HK_244" = "244",
  "HK_264" = "264", "HK_270" = "290",   # LR_270 = SR_290
  "PB_310" = "310", "PB_329" = "329", "PB_342" = "343",  # LR_342 = SR_343
  "PB_358" = "358", "PB_374" = "374", "PB_385" = "385"
)

# Ancestors to exclude from all analyses
SR_ANCESTORS <- c("HK104_ANC_1", "PB800_ANC__1", "anc1", "anc2")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Classify SNPs into 6 mutation classes
classify_snp <- function(ref, alt) {
  case_when(
    (ref == "G" & alt == "A") | (ref == "C" & alt == "T") ~ "GC_to_AT",
    (ref == "G" & alt == "C") | (ref == "C" & alt == "G") ~ "GC_to_CG",
    (ref == "G" & alt == "T") | (ref == "C" & alt == "A") ~ "GC_to_TA",
    (ref == "A" & alt == "G") | (ref == "T" & alt == "C") ~ "AT_to_GC",
    (ref == "A" & alt == "C") | (ref == "T" & alt == "G") ~ "AT_to_CG",
    (ref == "A" & alt == "T") | (ref == "T" & alt == "A") ~ "AT_to_TA",
    TRUE ~ NA_character_
  )
}

# Classify indels into 10 size bins
classify_indel <- function(ref, alt) {
  del_len <- nchar(ref) - nchar(alt)
  ins_len <- nchar(alt) - nchar(ref)
  case_when(
    del_len == 1                    ~ "del_1",
    del_len == 2                    ~ "del_2",
    del_len >= 3 & del_len <= 5     ~ "del_3_5",
    del_len >= 6 & del_len <= 10    ~ "del_6_10",
    del_len > 10                    ~ "del_10plus",
    ins_len == 1                    ~ "ins_1",
    ins_len == 2                    ~ "ins_2",
    ins_len >= 3 & ins_len <= 5     ~ "ins_3_5",
    ins_len >= 6 & ins_len <= 10    ~ "ins_6_10",
    ins_len > 10                    ~ "ins_10plus",
    TRUE ~ NA_character_
  )
}

# Compute per-site mutation rate
compute_rate <- function(count, callable_sites, generations = GENERATIONS) {
  ifelse(
    is.na(callable_sites) | callable_sites == 0 | is.na(count),
    NA_real_,
    count / (callable_sites * generations)
  )
}

# Summarise mutation counts and spectrum from a mutation data frame
# mut_df must have columns: Sample, REF, ALT, TYPE
# Returns one row per sample with counts for all 16 categories
summarise_mutations <- function(mut_df, sample_col = "Sample") {

  snp_classes  <- c("GC_to_AT","GC_to_CG","GC_to_TA","AT_to_GC","AT_to_CG","AT_to_TA")
  indel_bins   <- c("del_1","del_2","del_3_5","del_6_10","del_10plus",
                    "ins_1","ins_2","ins_3_5","ins_6_10","ins_10plus")

  mut_df <- mut_df %>%
    mutate(
      REF = as.character(REF),
      ALT = as.character(ALT),
      TYPE = tolower(as.character(TYPE)),
      snp_class   = ifelse(TYPE == "snp",   classify_snp(REF, ALT),   NA_character_),
      indel_class = ifelse(TYPE == "indel", classify_indel(REF, ALT), NA_character_)
    )

  # Per-sample totals
  totals <- mut_df %>%
    group_by(Sample = .data[[sample_col]]) %>%
    summarise(
      n_snp   = sum(TYPE == "snp",   na.rm = TRUE),
      n_indel = sum(TYPE == "indel", na.rm = TRUE),
      n_total = n(),
      .groups = "drop"
    )

  # SNP class counts (wide)
  snp_counts <- mut_df %>%
    filter(TYPE == "snp", !is.na(snp_class)) %>%
    group_by(Sample = .data[[sample_col]], snp_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = snp_class, values_from = n, values_fill = 0) %>%
    # ensure all 6 columns exist even if some classes absent
    { for (cl in snp_classes) if (!cl %in% names(.)) .[[cl]] <- 0L; . }

  # Indel bin counts (wide)
  indel_counts <- mut_df %>%
    filter(TYPE == "indel", !is.na(indel_class)) %>%
    group_by(Sample = .data[[sample_col]], indel_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = indel_class, values_from = n, values_fill = 0) %>%
    { for (bn in indel_bins) if (!bn %in% names(.)) .[[bn]] <- 0L; . }

  # Merge all
  result <- totals %>%
    left_join(snp_counts,   by = c("Sample" = "Sample")) %>%
    left_join(indel_counts, by = c("Sample" = "Sample"))

  # Fill any remaining NAs with 0 for count columns
  count_cols <- c(snp_classes, indel_bins)
  result[count_cols] <- lapply(result[count_cols], function(x) replace(x, is.na(x), 0))

  return(result)
}

# =============================================================================
# 1. LOAD METADATA
# =============================================================================
cat("Loading metadata...\n")

meta_sr_raw <- read_excel("C_briggsae_MetaData.xlsx", sheet = "Short Reads")
meta_lr_raw <- read_excel("C_briggsae_MetaData.xlsx", sheet = "Long Reads")

# SR metadata — one row per lane; take unique per line (concordant sites same for both lanes)
meta_sr <- meta_sr_raw %>%
  rename(
    Strain       = `Strain (HK104 or PB800)`,
    Line         = `Line (anc1, anc2, MA#)`,
    SR_callable_3x  = `# genome concordant sites (3x)`,
    SR_callable_10x = `# genome concordant sites (10x)`
  ) %>%
  filter(!Line %in% SR_ANCESTORS) %>%
  select(Strain, Line, SR_callable_3x, SR_callable_10x) %>%
  distinct(Line, .keep_all = TRUE) %>%
  mutate(Line = as.character(Line))

cat(sprintf("  SR metadata: %d lines\n", nrow(meta_sr)))

# LR metadata — one row per sample
meta_lr <- meta_lr_raw %>%
  rename(
    LR_sample       = `Strain (HK104 or PB800)`,
    LR_callable_3x  = `# genome covered 3X`,
    LR_callable_10x = `# genome covered 10X`
  ) %>%
  filter(!LR_sample %in% LR_EXCLUDE) %>%
  mutate(
    SR_Line = recode(LR_sample, !!!LR_TO_SR)
  ) %>%
  select(LR_sample, SR_Line, LR_callable_3x, LR_callable_10x)

cat(sprintf("  LR metadata: %d lines\n", nrow(meta_lr)))

# =============================================================================
# 2. LOAD SR MUTATIONS (all lines)
# =============================================================================
cat("Loading SR mutations...\n")

load_sr_sheet <- function(depth) {
  df <- read_excel("SR_final_mutation_list.xlsx", sheet = depth) %>%
    mutate(Sample = as.character(Sample)) %>%
    filter(!Sample %in% SR_ANCESTORS)
  cat(sprintf("  SR %s: %d mutations across %d lines\n",
              depth, nrow(df), n_distinct(df$Sample)))
  df
}

sr_3x  <- load_sr_sheet("3x")
sr_10x <- load_sr_sheet("10x")

# Summarise SR mutations
sr_summary_3x  <- summarise_mutations(sr_3x)  %>% rename_with(~paste0("SR_3x_",  .), -Sample)
sr_summary_10x <- summarise_mutations(sr_10x) %>% rename_with(~paste0("SR_10x_", .), -Sample)

# =============================================================================
# 3. LOAD SR-LR MUTATIONS (matched lines only)
# =============================================================================
cat("Loading SR-LR mutations...\n")

load_srlr_sheet <- function(depth) {
  df <- read_excel("SR_LR_final_mutation_list.xlsx", sheet = depth) %>%
    mutate(
      Sample = as.character(Sample),
      B1     = as.character(B1),
      B2     = as.character(B2),
      LR     = as.character(LR)
    ) %>%
    filter(!Sample %in% SR_ANCESTORS)
  cat(sprintf("  SR-LR %s: %d rows across %d lines\n",
              depth, nrow(df), n_distinct(df$Sample)))
  df
}

srlr_3x  <- load_srlr_sheet("3x")
srlr_10x <- load_srlr_sheet("10x")

# --- SR counts from the matched file (B1 or B2 not absent) ---
sr_matched_summary <- function(df, suffix) {
  df %>%
    filter(B1 != "absent" | B2 != "absent") %>%
    summarise_mutations() %>%
    rename_with(~paste0("SR_matched_", suffix, "_", .), -Sample)
}

# --- LR counts (LR not absent, regardless of SR) ---
lr_summary <- function(df, suffix) {
  df %>%
    filter(LR != "absent") %>%
    summarise_mutations() %>%
    rename_with(~paste0("LR_", suffix, "_", .), -Sample)
}

sr_matched_3x  <- sr_matched_summary(srlr_3x,  "3x")
sr_matched_10x <- sr_matched_summary(srlr_10x, "10x")
lr_sum_3x      <- lr_summary(srlr_3x,  "3x")
lr_sum_10x     <- lr_summary(srlr_10x, "10x")

# =============================================================================
# 4. BUILD MASTER TABLE
# =============================================================================
cat("Building master table...\n")

# Start from SR metadata (all lines)
master <- meta_sr %>%
  mutate(Generations = GENERATIONS) %>%
  # Join SR mutation summaries
  left_join(sr_summary_3x,  by = c("Line" = "Sample")) %>%
  left_join(sr_summary_10x, by = c("Line" = "Sample")) %>%
  # Join SR-matched summaries
  left_join(sr_matched_3x,  by = c("Line" = "Sample")) %>%
  left_join(sr_matched_10x, by = c("Line" = "Sample")) %>%
  # Join LR metadata (by SR line ID)
  left_join(meta_lr %>% select(SR_Line, LR_sample, LR_callable_3x, LR_callable_10x),
            by = c("Line" = "SR_Line")) %>%
  # Join LR mutation summaries
  left_join(lr_sum_3x,  by = c("Line" = "Sample")) %>%
  left_join(lr_sum_10x, by = c("Line" = "Sample"))

# =============================================================================
# 5. COMPUTE MUTATION RATES
# =============================================================================
cat("Computing mutation rates...\n")

master <- master %>%
  mutate(
    # SR rates (all lines)
    SR_SNP_rate_3x    = compute_rate(SR_3x_n_snp,   SR_callable_3x),
    SR_Indel_rate_3x  = compute_rate(SR_3x_n_indel, SR_callable_3x),
    SR_Total_rate_3x  = compute_rate(SR_3x_n_total, SR_callable_3x),
    SR_SNP_rate_10x   = compute_rate(SR_10x_n_snp,   SR_callable_10x),
    SR_Indel_rate_10x = compute_rate(SR_10x_n_indel, SR_callable_10x),
    SR_Total_rate_10x = compute_rate(SR_10x_n_total, SR_callable_10x),

    # SR matched rates (lines with LR data)
    SR_matched_SNP_rate_3x    = compute_rate(SR_matched_3x_n_snp,   SR_callable_3x),
    SR_matched_Indel_rate_3x  = compute_rate(SR_matched_3x_n_indel, SR_callable_3x),
    SR_matched_SNP_rate_10x   = compute_rate(SR_matched_10x_n_snp,  SR_callable_10x),
    SR_matched_Indel_rate_10x = compute_rate(SR_matched_10x_n_indel,SR_callable_10x),

    # LR rates (matched lines only)
    LR_SNP_rate_3x    = compute_rate(LR_3x_n_snp,   LR_callable_3x),
    LR_Indel_rate_3x  = compute_rate(LR_3x_n_indel, LR_callable_3x),
    LR_Total_rate_3x  = compute_rate(LR_3x_n_total, LR_callable_3x),
    LR_SNP_rate_10x   = compute_rate(LR_10x_n_snp,   LR_callable_10x),
    LR_Indel_rate_10x = compute_rate(LR_10x_n_indel, LR_callable_10x),
    LR_Total_rate_10x = compute_rate(LR_10x_n_total, LR_callable_10x),

    # Ts/Tv ratio (SR)
    SR_TsTv_3x  = (SR_3x_GC_to_AT  + SR_3x_AT_to_GC)  /
                  (SR_3x_GC_to_CG  + SR_3x_GC_to_TA  + SR_3x_AT_to_CG + SR_3x_AT_to_TA),
    SR_TsTv_10x = (SR_10x_GC_to_AT + SR_10x_AT_to_GC) /
                  (SR_10x_GC_to_CG + SR_10x_GC_to_TA + SR_10x_AT_to_CG + SR_10x_AT_to_TA),

    # Ts/Tv ratio (LR)
    LR_TsTv_3x  = (LR_3x_GC_to_AT  + LR_3x_AT_to_GC)  /
                  (LR_3x_GC_to_CG  + LR_3x_GC_to_TA  + LR_3x_AT_to_CG + LR_3x_AT_to_TA),
    LR_TsTv_10x = (LR_10x_GC_to_AT + LR_10x_AT_to_GC) /
                  (LR_10x_GC_to_CG + LR_10x_GC_to_TA + LR_10x_AT_to_CG + LR_10x_AT_to_TA),

    # Ins:Del ratio (SR)
    SR_InsDel_3x  = (SR_3x_ins_1  + SR_3x_ins_2  + SR_3x_ins_3_5  + SR_3x_ins_6_10  + SR_3x_ins_10plus) /
                    (SR_3x_del_1  + SR_3x_del_2  + SR_3x_del_3_5  + SR_3x_del_6_10  + SR_3x_del_10plus),
    SR_InsDel_10x = (SR_10x_ins_1 + SR_10x_ins_2 + SR_10x_ins_3_5 + SR_10x_ins_6_10 + SR_10x_ins_10plus) /
                    (SR_10x_del_1 + SR_10x_del_2 + SR_10x_del_3_5 + SR_10x_del_6_10 + SR_10x_del_10plus),

    # Ins:Del ratio (LR)
    LR_InsDel_3x  = (LR_3x_ins_1  + LR_3x_ins_2  + LR_3x_ins_3_5  + LR_3x_ins_6_10  + LR_3x_ins_10plus) /
                    (LR_3x_del_1  + LR_3x_del_2  + LR_3x_del_3_5  + LR_3x_del_6_10  + LR_3x_del_10plus),
    LR_InsDel_10x = (LR_10x_ins_1 + LR_10x_ins_2 + LR_10x_ins_3_5 + LR_10x_ins_6_10 + LR_10x_ins_10plus) /
                    (LR_10x_del_1 + LR_10x_del_2 + LR_10x_del_3_5 + LR_10x_del_6_10 + LR_10x_del_10plus)
  )

# =============================================================================
# 6. REORDER COLUMNS LOGICALLY
# =============================================================================
snp_classes  <- c("GC_to_AT","GC_to_CG","GC_to_TA","AT_to_GC","AT_to_CG","AT_to_TA")
indel_bins   <- c("del_1","del_2","del_3_5","del_6_10","del_10plus",
                  "ins_1","ins_2","ins_3_5","ins_6_10","ins_10plus")

id_cols    <- c("Line","Strain","Generations","LR_sample")
cov_cols   <- c("SR_callable_3x","SR_callable_10x","LR_callable_3x","LR_callable_10x")
rate_cols  <- c("SR_SNP_rate_3x","SR_Indel_rate_3x","SR_Total_rate_3x",
                "SR_SNP_rate_10x","SR_Indel_rate_10x","SR_Total_rate_10x",
                "SR_matched_SNP_rate_3x","SR_matched_Indel_rate_3x",
                "SR_matched_SNP_rate_10x","SR_matched_Indel_rate_10x",
                "LR_SNP_rate_3x","LR_Indel_rate_3x","LR_Total_rate_3x",
                "LR_SNP_rate_10x","LR_Indel_rate_10x","LR_Total_rate_10x")
diag_cols  <- c("SR_TsTv_3x","SR_TsTv_10x","LR_TsTv_3x","LR_TsTv_10x",
                "SR_InsDel_3x","SR_InsDel_10x","LR_InsDel_3x","LR_InsDel_10x")
count_cols <- c(
  paste0("SR_3x_n_snp"), paste0("SR_3x_n_indel"), paste0("SR_3x_n_total"),
  paste0("SR_10x_n_snp"),paste0("SR_10x_n_indel"),paste0("SR_10x_n_total"),
  paste0("SR_matched_3x_n_snp"), paste0("SR_matched_3x_n_indel"),
  paste0("SR_matched_10x_n_snp"),paste0("SR_matched_10x_n_indel"),
  paste0("LR_3x_n_snp"), paste0("LR_3x_n_indel"), paste0("LR_3x_n_total"),
  paste0("LR_10x_n_snp"),paste0("LR_10x_n_indel"),paste0("LR_10x_n_total")
)
spec_cols  <- c(
  paste0("SR_3x_",  snp_classes), paste0("SR_3x_",  indel_bins),
  paste0("SR_10x_", snp_classes), paste0("SR_10x_", indel_bins),
  paste0("SR_matched_3x_",  snp_classes), paste0("SR_matched_3x_",  indel_bins),
  paste0("SR_matched_10x_", snp_classes), paste0("SR_matched_10x_", indel_bins),
  paste0("LR_3x_",  snp_classes), paste0("LR_3x_",  indel_bins),
  paste0("LR_10x_", snp_classes), paste0("LR_10x_", indel_bins)
)

# Select only columns that exist
all_cols <- c(id_cols, cov_cols, rate_cols, diag_cols, count_cols, spec_cols)
all_cols <- all_cols[all_cols %in% names(master)]

master <- master %>% select(all_of(all_cols))

# =============================================================================
# 7. SAVE
# =============================================================================
write.csv(master, "master_table.csv", row.names = FALSE)
cat(sprintf("\nMaster table saved: %d lines x %d columns\n", nrow(master), ncol(master)))

# Quick summary
cat("\n--- SR-only summary (3x) ---\n")
master %>%
  group_by(Strain) %>%
  summarise(
    n_lines        = n(),
    mean_SNP_rate  = mean(SR_SNP_rate_3x,   na.rm = TRUE),
    mean_Indel_rate= mean(SR_Indel_rate_3x, na.rm = TRUE),
    mean_TsTv      = mean(SR_TsTv_3x,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

cat("\n--- SR vs LR (3x, matched lines only) ---\n")
master %>%
  filter(!is.na(LR_sample)) %>%
  group_by(Strain) %>%
  summarise(
    n_lines            = n(),
    SR_mean_SNP_rate   = mean(SR_matched_SNP_rate_3x, na.rm = TRUE),
    LR_mean_SNP_rate   = mean(LR_SNP_rate_3x,         na.rm = TRUE),
    SR_mean_TsTv       = mean(SR_TsTv_3x,              na.rm = TRUE),
    LR_mean_TsTv       = mean(LR_TsTv_3x,              na.rm = TRUE),
    SR_mean_InsDel     = mean(SR_InsDel_3x,            na.rm = TRUE),
    LR_mean_InsDel     = mean(LR_InsDel_3x,            na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

cat("\nDone.\n")
