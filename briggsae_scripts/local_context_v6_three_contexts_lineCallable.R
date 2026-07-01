# =============================================================================
# Local Sequence Context (Trinucleotide) Analysis -- C. briggsae MA Lines
# v6: three context orientations + per-line callable denominator
#
# Major v6 change:
#   Rates now use each line's own callable genome size from master_table:
#
#       line_callable_fraction_3x = SR_callable_3x / GENOME_SIZE
#
#       rate_line_motif =
#           n_mut_line_motif /
#           (genome_motif_count * line_callable_fraction_3x * line_generations)
#
#   This replaces the older global denominator:
#
#       genome_motif_count * GENERATIONS * CALLABLE_FRAC
#
# Important limitation:
#   This is still line-specific GLOBAL callability, not motif-specific
#   callable counts. The ideal denominator would be motif-specific callable
#   counts per line, but master_table only gives total callable sites per line.
#
# Context orientations:
#   1) Centered   : 5'-xYz-3'   (Rajaei-style; Y = mutable middle base)
#   2) 3-prime    : 5'-xxY-3'   (Saxena-style; Y = mutable 3' base)
#   3) 5-prime    : 5'-Yxx-3'   (complementary check; Y = mutable 5' base)
#
# Required inputs:
#   SR_final_mutation_list_mono_classified.xlsx
#   briggsae_genome_trinucleotide_counts.csv
#   master_table.csv OR master_table(1).csv
#   mutations_snp_3bp_center_xYz.fasta
#   mutations_snp_3bp_3prime_xxY.fasta
#   mutations_snp_3bp_5prime_Yxx.fasta
# =============================================================================

library(rstudioapi)
library(Biostrings)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd(dirname(getActiveDocumentContext()$path))

# =============================================================================
# CONSTANTS
# =============================================================================
GENOME_SIZE <- 106196309

# Kept only as a reminder/fallback reference. Not used in v6 rate calculations.
OLD_GLOBAL_CALLABLE_FRAC <- 0.91

SUB_COLORS <- c(
  "C>A" = "#03BCEE", "C>G" = "#010101", "C>T" = "#E32926",
  "T>A" = "#CAC9C9", "T>C" = "#A1CE63", "T>G" = "#ECC7C4"
)

NUC_COLORS <- c(A = "#00BB00", C = "#0000BB", G = "#E47E00", T = "#BB0000")

context_configs <- tibble::tribble(
  ~context_id,     ~fasta_file,                           ~focal_idx, ~context_label,                  ~focal_name,    ~outfile_stub,
  "center_xYz",    "mutations_snp_3bp_center_xYz.fasta",   2L,         "5'-xYz-3'; Y = middle base",   "Middle base",  "center_xYz",
  "3prime_xxY",    "mutations_snp_3bp_3prime_xxY.fasta",   3L,         "5'-xxY-3'; Y = 3' base",       "3' base",      "3prime_xxY",
  "5prime_Yxx",    "mutations_snp_3bp_5prime_Yxx.fasta",   1L,         "5'-Yxx-3'; Y = 5' base",       "5' base",      "5prime_Yxx"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
save_fig <- function(p, name, w, h) {
  ggsave(paste0(name, ".pdf"), p, width = w, height = h)
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s.pdf/.png\n", name))
}

theme_tri <- function() {
  theme_bw(base_size = 20) +
    theme(
      axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        size = 10, family = "mono", colour = "black", face = "bold"
      ),
      axis.text.y      = element_text(size = 18, colour = "black"),
      axis.title       = element_text(size = 22),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(size = 20, face = "bold"),
      legend.position  = "none",
      panel.spacing    = unit(0.1, "lines"),
      plot.title       = element_blank(),
      plot.subtitle    = element_text(size = 22, face = "bold")
    )
}

theme_tri64 <- function() {
  theme_bw(base_size = 20) +
    theme(
      axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        size = 10, family = "mono", colour = "black", face = "bold"
      ),
      axis.text.y      = element_text(size = 18, colour = "black"),
      axis.title       = element_text(size = 22),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(size = 20, face = "bold"),
      legend.position  = "top",
      legend.title     = element_text(size = 20),
      legend.text      = element_text(size = 20),
      panel.spacing    = unit(0.12, "lines"),
      plot.title       = element_blank(),
      plot.subtitle    = element_text(size = 24, face = "bold")
    )
}

ylab_rate <- expression("Mean " * mu ~ "(" * 10^{-9} ~ "per site per gen)")

make_motif_levels <- function(focal_idx) {
  b <- c("A", "C", "G", "T")

  if (focal_idx == 2L) {
    expand.grid(x = b, Y = b, z = b, stringsAsFactors = FALSE) %>%
      arrange(Y, x, z) %>%
      mutate(
        context_64 = paste0(x, Y, z),
        label_64   = paste0(tolower(x), Y, tolower(z)),
        focal_base = Y
      )
  } else if (focal_idx == 3L) {
    expand.grid(x1 = b, x2 = b, Y = b, stringsAsFactors = FALSE) %>%
      arrange(Y, x1, x2) %>%
      mutate(
        context_64 = paste0(x1, x2, Y),
        label_64   = paste0(tolower(x1), tolower(x2), Y),
        focal_base = Y
      )
  } else if (focal_idx == 1L) {
    expand.grid(Y = b, x1 = b, x2 = b, stringsAsFactors = FALSE) %>%
      arrange(Y, x1, x2) %>%
      mutate(
        context_64 = paste0(Y, x1, x2),
        label_64   = paste0(Y, tolower(x1), tolower(x2)),
        focal_base = Y
      )
  } else {
    stop("focal_idx must be 1, 2, or 3")
  }
}

parse_fasta_context <- function(fasta_file, focal_idx) {
  cat(sprintf("Parsing FASTA: %s\n", fasta_file))

  fasta <- readDNAStringSet(fasta_file)

  fasta_df <- data.frame(
    fasta_id = names(fasta),
    context  = toupper(as.character(fasta)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      mut_id = sub("::.*", "", fasta_id)
    ) %>%
    filter(nchar(context) == 3, !grepl("N", context))

  # Parse header using a right-anchored expression:
  # CHROM_POS_REF_ALT_Sample
  # This is safer than splitting from the left if chromosome names contain underscores.
  m <- regexec("^(.+)_(\\d+)_([ACGTacgt])_([ACGTacgt])_(.+)$", fasta_df$mut_id)
  parsed <- regmatches(fasta_df$mut_id, m)

  parsed_mat <- do.call(
    rbind,
    lapply(parsed, function(x) {
      if (length(x) == 6) x else rep(NA_character_, 6)
    })
  )

  colnames(parsed_mat) <- c("full", "CHROM", "POS", "REF", "ALT", "Sample")

  fasta_df <- fasta_df %>%
    bind_cols(as.data.frame(parsed_mat[, -1], stringsAsFactors = FALSE)) %>%
    mutate(
      POS = as.integer(POS),
      REF = toupper(REF),
      ALT = toupper(ALT),
      sample_num = suppressWarnings(as.integer(sub("_.*", "", as.character(Sample)))),
      focal_base = substr(context, focal_idx, focal_idx),
      focal_matches_ref = focal_base == REF
    )

  n_bad <- sum(!fasta_df$focal_matches_ref, na.rm = TRUE)
  n_non_numeric <- sum(is.na(fasta_df$sample_num))

  cat(sprintf("  Parsed %d usable 3-bp contexts\n", nrow(fasta_df)))
  cat(sprintf("  Focal base != REF: %d\n", n_bad))
  cat(sprintf("  Non-numeric sample IDs in FASTA: %d\n", n_non_numeric))

  fasta_df
}

add_line_coverage <- function(line_strain, master_cov) {
  out <- line_strain %>%
    left_join(master_cov, by = "sample_num")

  missing_cov <- out %>%
    filter(
      is.na(line_callable_sites_3x) |
      is.na(line_callable_fraction_3x) |
      is.na(line_generations)
    )

  if (nrow(missing_cov) > 0) {
    print(missing_cov)
    stop("Some lines are missing SR_callable_3x or Generations in master_table. Fix master_table or Sample/Line matching.")
  }

  out
}

compute_64_rates <- function(muts_base, fasta_df, tri_genome_64, master_cov,
                             focal_idx, context_id, context_label, focal_name) {
  motif_levels <- make_motif_levels(focal_idx)

  muts_local <- muts_base %>%
    left_join(fasta_df %>% select(mut_id, context),
              by = "mut_id") %>%
    filter(!is.na(context))

  line_strain <- muts_local %>%
    distinct(Sample, sample_num, Strain) %>%
    add_line_coverage(master_cov)

  full_grid_64 <- line_strain %>%
    cross_join(tri_genome_64 %>% select(context_64, genome_count_64)) %>%
    left_join(
      muts_local %>% count(Sample, context_64 = context, name = "n_mut_64"),
      by = c("Sample", "context_64")
    ) %>%
    replace_na(list(n_mut_64 = 0)) %>%
    mutate(
      # v6 denominator: line-specific callable fraction and line-specific generations
      effective_callable_motif_count = genome_count_64 * line_callable_fraction_3x,
      rate_64 = n_mut_64 / (effective_callable_motif_count * line_generations)
    )

  mut_counts_64 <- full_grid_64 %>%
    group_by(Strain, context_64) %>%
    summarise(
      n_lines      = n(),
      n_mut        = sum(n_mut_64),
      mean_rate    = mean(rate_64),
      sd_rate      = sd(rate_64),
      sem_rate     = sd(rate_64) / sqrt(n()),
      mean_callable_fraction_3x = mean(line_callable_fraction_3x),
      min_callable_fraction_3x  = min(line_callable_fraction_3x),
      max_callable_fraction_3x  = max(line_callable_fraction_3x),
      .groups      = "drop"
    ) %>%
    left_join(tri_genome_64, by = "context_64") %>%
    left_join(motif_levels, by = "context_64") %>%
    mutate(
      label_64   = factor(label_64, levels = motif_levels$label_64),
      context_id = context_id,
      context_label = context_label,
      focal_name = focal_name
    )

  list(muts_local = muts_local, mut_counts_64 = mut_counts_64, motif_levels = motif_levels)
}

plot_64motif <- function(mut_counts_64, strain_name, context_label, focal_name) {
  d <- mut_counts_64 %>%
    filter(Strain == strain_name) %>%
    arrange(focal_base, label_64)

  ggplot(d, aes(x = label_64, y = mean_rate * 1e9, fill = focal_base)) +
    geom_col(width = 0.8, color = NA) +
    geom_errorbar(
      aes(ymin = pmax(0, (mean_rate - sem_rate) * 1e9),
          ymax = (mean_rate + sem_rate) * 1e9),
      width = 0.35, linewidth = 0.35, colour = "grey30"
    ) +
    scale_fill_manual(values = NUC_COLORS, name = focal_name) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_grid(. ~ focal_base, scales = "free_x", space = "free_x") +
    labs(
      x = paste0("64 three-bp motifs (", context_label, ")"),
      y = ylab_rate,
      subtitle = strain_name
    ) +
    theme_tri64()
}

# =============================================================================
# STEP 1: LOAD MUTATION LIST
# =============================================================================
cat("Loading mutation list...\n")

muts_base <- read_excel("SR_final_mutation_list_mono_classified.xlsx", sheet = "3x") %>%
  filter(TYPE == "snp") %>%
  mutate(
    REF = toupper(as.character(REF)),
    ALT = toupper(as.character(ALT)),
    sample_num = suppressWarnings(as.integer(sub("_.*", "", as.character(Sample)))),
    Strain = case_when(
      sample_num <= 298 ~ "HK104",
      sample_num >= 300 ~ "PB800",
      TRUE ~ NA_character_
    ),
    mut_id = paste(CHROM, POS, REF, ALT, Sample, sep = "_")
  ) %>%
  filter(!is.na(Strain))

cat(sprintf("  %d SNV mutations loaded (%d HK104, %d PB800)\n",
            nrow(muts_base), sum(muts_base$Strain == "HK104"), sum(muts_base$Strain == "PB800")))

# =============================================================================
# STEP 2: LOAD MASTER TABLE WITH PER-LINE CALLABILITY
# =============================================================================
master_candidates <- c("master_table.csv", "master_table(1).csv")
master_file <- master_candidates[file.exists(master_candidates)][1]

if (is.na(master_file)) {
  stop("Could not find master_table.csv or master_table(1).csv in the working directory.")
}

cat(sprintf("Loading per-line callable information from: %s\n", master_file))

master_cov <- read.csv(master_file, stringsAsFactors = FALSE) %>%
  transmute(
    sample_num = as.integer(Line),
    master_strain = as.character(Strain),
    line_generations = as.numeric(Generations),
    line_callable_sites_3x = as.numeric(SR_callable_3x),
    line_callable_sites_10x = as.numeric(SR_callable_10x),
    line_callable_fraction_3x = line_callable_sites_3x / GENOME_SIZE,
    line_callable_fraction_10x = line_callable_sites_10x / GENOME_SIZE
  )

cat("Per-line SR callable fraction summary, 3x:\n")
print(summary(master_cov$line_callable_fraction_3x))

write.csv(master_cov, "line_callable_denominators_from_master_table.csv", row.names = FALSE)
cat("Saved: line_callable_denominators_from_master_table.csv\n")

# =============================================================================
# STEP 3: GENOME TRINUCLEOTIDE COUNTS
# =============================================================================
cat("Loading genome trinucleotide counts...\n")

comp <- c(A = "T", T = "A", C = "G", G = "C",
          a = "T", t = "A", c = "G", g = "C")

tri_genome <- read.csv("briggsae_genome_trinucleotide_counts.csv",
                       stringsAsFactors = FALSE) %>%
  mutate(
    trinucleotide = toupper(trinucleotide),
    mid     = substr(trinucleotide, 2, 2),
    ctx_pyr = ifelse(mid %in% c("A", "G"),
                     paste0(comp[substr(trinucleotide, 3, 3)],
                            comp[substr(trinucleotide, 2, 2)],
                            comp[substr(trinucleotide, 1, 1)]),
                     trinucleotide)
  ) %>%
  group_by(context_pyr = ctx_pyr) %>%
  summarise(genome_count = sum(count), .groups = "drop")

tri_genome_64 <- read.csv("briggsae_genome_trinucleotide_counts.csv",
                          stringsAsFactors = FALSE) %>%
  mutate(trinucleotide = toupper(trinucleotide)) %>%
  filter(nchar(trinucleotide) == 3, grepl("^[ACGT]{3}$", trinucleotide)) %>%
  group_by(context_64 = trinucleotide) %>%
  summarise(genome_count_64 = sum(count), .groups = "drop")

# =============================================================================
# STEP 4: ORIGINAL 96-CHANNEL ANALYSIS FOR CENTERED CONTEXT ONLY
# =============================================================================
cat("\n=== Centered 96-channel analysis with line-specific callable denominator ===\n")

center_cfg <- context_configs %>% filter(context_id == "center_xYz")
center_fasta <- parse_fasta_context(center_cfg$fasta_file[[1]], center_cfg$focal_idx[[1]])

muts_center <- muts_base %>%
  left_join(center_fasta %>% select(mut_id, context),
            by = "mut_id")

n_missing <- sum(is.na(muts_center$context))
cat(sprintf("  Centered context joined: %d / %d (%d missing)\n",
            nrow(muts_center) - n_missing, nrow(muts_center), n_missing))

muts_center <- muts_center %>% filter(!is.na(context))

muts_96 <- muts_center %>%
  filter(!is.na(REF), !is.na(ALT), !is.na(context),
         nchar(context) == 3, !grepl("N", context)) %>%
  mutate(
    REF = toupper(REF), ALT = toupper(ALT), ctx = toupper(context),
    is_purine = REF %in% c("A", "G"),
    ctx_pyr = ifelse(is_purine,
                     paste0(comp[substr(ctx, 3, 3)], comp[substr(ctx, 2, 2)], comp[substr(ctx, 1, 1)]),
                     ctx),
    alt_pyr  = ifelse(is_purine, comp[ALT], ALT),
    ref_pyr  = ifelse(is_purine, comp[REF], REF),
    sub_type = paste0(ref_pyr, ">", alt_pyr),
    context_pyr = ctx_pyr
  ) %>%
  select(-is_purine, -ctx, -ctx_pyr, -alt_pyr, -ref_pyr)

line_strain <- muts_96 %>%
  distinct(Sample, sample_num, Strain) %>%
  add_line_coverage(master_cov)

ctx_sub_combos <- muts_96 %>% distinct(context_pyr, sub_type)

full_grid <- line_strain %>%
  cross_join(ctx_sub_combos) %>%
  left_join(
    muts_96 %>% count(Sample, context_pyr, sub_type, name = "n_mut"),
    by = c("Sample", "context_pyr", "sub_type")
  ) %>%
  replace_na(list(n_mut = 0)) %>%
  left_join(tri_genome, by = "context_pyr") %>%
  mutate(
    # v6 denominator: line-specific callable fraction and line-specific generations
    effective_callable_context_count = genome_count * line_callable_fraction_3x,
    rate = n_mut / (effective_callable_context_count * line_generations)
  )

mut_counts <- full_grid %>%
  group_by(Strain, context_pyr, sub_type) %>%
  summarise(
    n_lines   = n(),
    n_mut     = sum(n_mut),
    mean_rate = mean(rate),
    sd_rate   = sd(rate),
    sem_rate  = sd(rate) / sqrt(n()),
    mean_callable_fraction_3x = mean(line_callable_fraction_3x),
    min_callable_fraction_3x  = min(line_callable_fraction_3x),
    max_callable_fraction_3x  = max(line_callable_fraction_3x),
    .groups   = "drop"
  ) %>%
  left_join(tri_genome %>% select(context_pyr, genome_count),
            by = "context_pyr") %>%
  mutate(
    x     = substr(context_pyr, 1, 1),
    Y     = substr(context_pyr, 2, 2),
    z     = substr(context_pyr, 3, 3),
    label = paste0(x, "[", sub_type, "]", z)
  )

write.csv(mut_counts, "trinucleotide_rates_lineCallable.csv", row.names = FALSE)
cat("Saved: trinucleotide_rates_lineCallable.csv\n")

plot_trinucleotide <- function(strain_name) {
  d <- mut_counts %>%
    filter(Strain == strain_name) %>%
    arrange(sub_type, x, z) %>%
    mutate(label = factor(label, levels = unique(label)))

  ggplot(d, aes(x = label, y = mean_rate * 1e9, fill = sub_type)) +
    geom_col(width = 0.8, color = NA) +
    geom_errorbar(
      aes(ymin = pmax(0, (mean_rate - sem_rate) * 1e9),
          ymax = (mean_rate + sem_rate) * 1e9),
      width = 0.4, linewidth = 0.3, colour = "grey30"
    ) +
    scale_fill_manual(values = SUB_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_grid(. ~ sub_type, scales = "free_x", space = "free_x") +
    labs(x = "Trinucleotide context", y = ylab_rate, subtitle = strain_name) +
    theme_tri()
}

p_hk <- plot_trinucleotide("HK104")
p_pb <- plot_trinucleotide("PB800")

save_fig(p_hk, "Fig_trinucleotide_HK104_lineCallable", 14, 4.5)
save_fig(p_pb, "Fig_trinucleotide_PB800_lineCallable", 14, 4.5)

fig_combined <- p_hk / p_pb
save_fig(fig_combined, "Fig_trinucleotide_combined_lineCallable", 14, 9)

# =============================================================================
# STEP 5: 64-MOTIF ANALYSIS FOR ALL THREE CONTEXT TYPES
# =============================================================================
cat("\n=== 64-motif exact-context analysis for all three context types ===\n")
cat("Denominator: genome motif count x line-specific SR_callable_3x fraction x line-specific generations\n")

all_64 <- list()
validation_rows <- list()

for (i in seq_len(nrow(context_configs))) {
  cfg <- context_configs[i, ]
  cat(sprintf("\n--- Context: %s ---\n", cfg$context_label[[1]]))

  fasta_df <- parse_fasta_context(cfg$fasta_file[[1]], cfg$focal_idx[[1]])

  validation_rows[[i]] <- data.frame(
    context_id = cfg$context_id[[1]],
    fasta_file = cfg$fasta_file[[1]],
    focal_idx  = cfg$focal_idx[[1]],
    n_records_usable_3bp_noN  = nrow(fasta_df),
    focal_base_not_ref = sum(!fasta_df$focal_matches_ref, na.rm = TRUE),
    non_numeric_sample_ids_in_fasta = sum(is.na(fasta_df$sample_num)),
    stringsAsFactors = FALSE
  )

  res <- compute_64_rates(
    muts_base     = muts_base,
    fasta_df      = fasta_df,
    tri_genome_64 = tri_genome_64,
    master_cov    = master_cov,
    focal_idx     = cfg$focal_idx[[1]],
    context_id    = cfg$context_id[[1]],
    context_label = cfg$context_label[[1]],
    focal_name    = cfg$focal_name[[1]]
  )

  mut_counts_64 <- res$mut_counts_64
  all_64[[cfg$context_id[[1]]]] <- mut_counts_64

  write.csv(
    mut_counts_64,
    paste0("trinucleotide_64motif_rates_", cfg$outfile_stub[[1]], "_lineCallable.csv"),
    row.names = FALSE
  )
  cat(sprintf("Saved: trinucleotide_64motif_rates_%s_lineCallable.csv\n", cfg$outfile_stub[[1]]))

  p64_hk <- plot_64motif(mut_counts_64, "HK104", cfg$context_label[[1]], cfg$focal_name[[1]])
  p64_pb <- plot_64motif(mut_counts_64, "PB800", cfg$context_label[[1]], cfg$focal_name[[1]])

  fig_64_combined <- p64_hk / p64_pb +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 18)
    )

  save_fig(
    fig_64_combined,
    paste0("Fig_trinucleotide_64motif_combined_", cfg$outfile_stub[[1]], "_lineCallable"),
    16, 9.5
  )

  save_fig(
    p64_hk,
    paste0("Fig_trinucleotide_64motif_HK104_", cfg$outfile_stub[[1]], "_lineCallable"),
    16, 4.75
  )

  save_fig(
    p64_pb,
    paste0("Fig_trinucleotide_64motif_PB800_", cfg$outfile_stub[[1]], "_lineCallable"),
    16, 4.75
  )
}

validation_df <- bind_rows(validation_rows)
write.csv(validation_df, "trinucleotide_64motif_context_validation_lineCallable.csv", row.names = FALSE)
cat("Saved: trinucleotide_64motif_context_validation_lineCallable.csv\n")

all_64_df <- bind_rows(all_64)
write.csv(all_64_df, "trinucleotide_64motif_rates_all_contexts_lineCallable.csv", row.names = FALSE)
cat("Saved: trinucleotide_64motif_rates_all_contexts_lineCallable.csv\n")

# =============================================================================
# STEP 6: QUICK SUMMARY
# =============================================================================
cat("\n=== TOP 8 64-MOTIF CONTEXTS PER STRAIN / CONTEXT TYPE ===\n")
top_64 <- all_64_df %>%
  filter(!is.na(mean_rate), mean_rate > 0) %>%
  group_by(context_id, context_label, Strain) %>%
  slice_max(mean_rate, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(mean_1e9 = mean_rate * 1e9,
         sem_1e9 = sem_rate * 1e9) %>%
  select(context_id, context_label, Strain, label_64, context_64, n_mut, mean_1e9, sem_1e9)

print(top_64)
write.csv(top_64, "trinucleotide_64motif_top8_all_contexts_lineCallable.csv", row.names = FALSE)
cat("Saved: trinucleotide_64motif_top8_all_contexts_lineCallable.csv\n")

cat("\nDone.\n")
cat("v6 denominator used: genome motif count x per-line SR_callable_3x fraction x per-line Generations.\n")
cat("Global CALLABLE_FRAC = 0.91 was not used in the v6 rate calculations.\n")
