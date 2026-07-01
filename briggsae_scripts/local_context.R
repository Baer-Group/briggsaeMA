# =============================================================================
# Local Sequence Context (Trinucleotide) Analysis -- C. briggsae MA Lines
# Equivalent of Rajaei et al. 2021 Figure 3
# =============================================================================
#
# WORKFLOW OVERVIEW:
#   Bash (run on HiPerGator first):
#     1. Convert mutation list to BED with 3-bp window
#     2. bedtools getfasta to extract trinucleotide context per mutation
#     3. Python to count trinucleotides in reference genome (normalization)
#
#   R (this script):
#     1. Parse fasta contexts + assign mutation class
#     2. Collapse to pyrimidine convention (32 unique contexts)
#     3. Compute per-line rates, then mean +/- SEM across lines
#     4. Plot 64-bar figure per strain (Rajaei Fig. 3 equivalent)
#     5. Correlation between HK104 and PB800 context rates
#
# Inputs:
#   SR_final_mutation_list_mono_classified.xlsx  (mutation list, 3x sheet)
#   mutations_snp_3bp.fasta     (output of bedtools getfasta -- bash step 2)
#   briggsae_genome_trinucleotide_counts.csv  (output of python -- bash step 3)
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
GENOME_SIZE   <- 106196309
CALLABLE_FRAC <- 0.91
GENERATIONS   <- 238
GC_FRAC       <- 0.37
AT_FRAC       <- 0.63
N_HK104       <- 32
N_PB800       <- 49

SUB_COLORS <- c(
  "C>A" = "#03BCEE", "C>G" = "#010101", "C>T" = "#E32926",
  "T>A" = "#CAC9C9", "T>C" = "#A1CE63", "T>G" = "#ECC7C4"
)

# =============================================================================
# BASH COMMANDS (run on HiPerGator before this script)
# =============================================================================
cat("
=== BASH COMMANDS TO RUN ON HIPERGATOR FIRST ===

# 1. Create BED with 3-bp window around each SNV
awk 'NR>1 && $8==\"snp\" {
    print $1\"\\t\"($2-2)\"\\t\"($2+1)\"\\t\"$1\"_\"$2\"_\"$4\"_\"$5\"_\"$6
}' SR_mutations_3x.tsv > mutations_snp_3bp.bed

# 2. Extract 3-bp fasta context
bedtools getfasta \\
    -fi 20250626_c_briggsae_Feb2020.genome.fa \\
    -bed mutations_snp_3bp.bed \\
    -fo mutations_snp_3bp.fasta \\
    -name

# 3. Count all trinucleotides in reference genome
#    Run briggsae_count_trinucleotides.py

===============================================
")

# =============================================================================
# STEP 1: LOAD MUTATION LIST
# =============================================================================
cat("Loading mutation list...\n")

muts <- read_excel("SR_final_mutation_list_mono_classified.xlsx", sheet = "3x") %>%
  filter(TYPE == "snp") %>%
  mutate(
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
            nrow(muts), sum(muts$Strain == "HK104"), sum(muts$Strain == "PB800")))

# =============================================================================
# STEP 2: PARSE FASTA CONTEXTS
# Header format: CHROM_POS_REF_ALT_Sample::start-end
# rest = POS_REF_ALT_Sample after removing CHROM_
# =============================================================================
cat("Parsing trinucleotide contexts from fasta...\n")

fasta <- readDNAStringSet("mutations_snp_3bp.fasta")
fasta_df <- data.frame(
  fasta_id = names(fasta),
  context  = toupper(as.character(fasta)),
  stringsAsFactors = FALSE
) %>%
  filter(nchar(context) == 3, !grepl("N", context))

fasta_df <- fasta_df %>%
  mutate(
    header_clean = sub("::.*", "", fasta_id),
    CHROM = sub("_.*", "", header_clean),
    rest  = sub("^[^_]+_", "", header_clean),
    POS   = as.integer(sub("_.*", "", rest)),
    REF   = substr(sub("^[^_]+_", "", rest), 1, 1),
    ALT   = substr(sub("^[^_]+_[^_]+_", "", rest), 1, 1)
  )

muts <- muts %>%
  left_join(fasta_df %>% select(CHROM, POS, REF, ALT, context),
            by = c("CHROM", "POS", "REF", "ALT"))

n_missing <- sum(is.na(muts$context))
cat(sprintf("  Context: %d / %d (%d missing)\n",
            nrow(muts) - n_missing, nrow(muts), n_missing))
muts <- muts %>% filter(!is.na(context))

# =============================================================================
# STEP 3: PYRIMIDINE CONVENTION
# If REF is A or G: reverse-complement context, REF, and ALT
# Result: REF is always C or T
# =============================================================================
comp <- c(A = "T", T = "A", C = "G", G = "C",
          a = "T", t = "A", c = "G", g = "C")

muts <- muts %>%
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

# =============================================================================
# STEP 4: GENOME TRINUCLEOTIDE COUNTS
# =============================================================================
cat("Loading genome trinucleotide counts...\n")

tri_genome <- read.csv("briggsae_genome_trinucleotide_counts.csv",
                       stringsAsFactors = FALSE) %>%
  mutate(
    mid     = substr(trinucleotide, 2, 2),
    ctx_pyr = ifelse(mid %in% c("A", "G"),
                     paste0(comp[substr(trinucleotide, 3, 3)],
                            comp[substr(trinucleotide, 2, 2)],
                            comp[substr(trinucleotide, 1, 1)]),
                     trinucleotide)
  ) %>%
  group_by(context_pyr = ctx_pyr) %>%
  summarise(genome_count = sum(count), .groups = "drop") %>%
  mutate(callable_count = genome_count * CALLABLE_FRAC)

cat(sprintf("  %d unique pyrimidine contexts\n", nrow(tri_genome)))

# =============================================================================
# STEP 5: PER-LINE RATES -> MEAN + SEM
# mu_i = n_mut_line_i / (callable_count * GENERATIONS)
# Then mean(mu_i) and SEM = sd(mu_i) / sqrt(n_lines) per strain
# =============================================================================
cat("Computing per-line rates with SEM...\n")

line_strain    <- muts %>% distinct(Sample, Strain)
ctx_sub_combos <- muts %>% distinct(context_pyr, sub_type)

# Full grid: every line x every (context, sub_type), zeros where no mutations
full_grid <- line_strain %>%
  cross_join(ctx_sub_combos) %>%
  left_join(
    muts %>% count(Sample, context_pyr, sub_type, name = "n_mut"),
    by = c("Sample", "context_pyr", "sub_type")
  ) %>%
  replace_na(list(n_mut = 0)) %>%
  left_join(tri_genome, by = "context_pyr") %>%
  mutate(rate = n_mut / (callable_count * GENERATIONS))

# Summarise across lines within each strain
mut_counts <- full_grid %>%
  group_by(Strain, context_pyr, sub_type) %>%
  summarise(
    n_lines   = n(),
    n_mut     = sum(n_mut),
    mean_rate = mean(rate),
    sd_rate   = sd(rate),
    sem_rate  = sd(rate) / sqrt(n()),
    .groups   = "drop"
  ) %>%
  left_join(tri_genome %>% select(context_pyr, genome_count, callable_count),
            by = "context_pyr") %>%
  mutate(
    x     = substr(context_pyr, 1, 1),
    Y     = substr(context_pyr, 2, 2),
    z     = substr(context_pyr, 3, 3),
    label = paste0(x, "[", sub_type, "]", z)
  )

write.csv(mut_counts, "trinucleotide_rates.csv", row.names = FALSE)
cat("Saved: trinucleotide_rates.csv\n")

# =============================================================================
# STEP 6: FIGURES
# =============================================================================

# -- Theme: larger text, bold black readable x-axis labels --
theme_tri <- function() {
  theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        size = 8, family = "mono", colour = "black", face = "bold"
      ),
      axis.text.y      = element_text(size = 12, colour = "black"),
      axis.title       = element_text(size = 14),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(size = 13, face = "bold"),
      legend.position  = "none",
      panel.spacing    = unit(0.1, "lines"),
      plot.title       = element_blank(),
      plot.subtitle    = element_text(size = 14, face = "bold")
    )
}

save_fig <- function(p, name, w, h) {
  ggsave(paste0(name, ".pdf"), p, width = w, height = h)
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s\n", name))
}

# Y-axis label: pure plotmath, no Unicode
ylab_rate <- expression("Mean " * mu ~ "(" * 10^{-9} ~ "per site per gen)")

# -- Fig 1: 64-bar trinucleotide context with SEM error bars --
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

save_fig(p_hk, "Fig_trinucleotide_HK104", 14, 4.5)
save_fig(p_pb, "Fig_trinucleotide_PB800", 14, 4.5)

fig_combined <- p_hk / p_pb
save_fig(fig_combined, "Fig_trinucleotide_combined", 14, 9)

# -- Fig 2: Correlation HK104 vs PB800 --
corr_data <- mut_counts %>%
  select(context_pyr, sub_type, Strain, mean_rate) %>%
  pivot_wider(names_from = Strain, values_from = mean_rate, values_fill = 0)

r_val <- cor(corr_data$HK104, corr_data$PB800, use = "complete.obs")
cat(sprintf("\nCorrelation HK104 vs PB800: r = %.4f\n", r_val))

fig_corr <- ggplot(corr_data,
                   aes(x = HK104 * 1e9, y = PB800 * 1e9, color = sub_type)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = SUB_COLORS, name = "Substitution") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %.3f", r_val),
           hjust = -0.2, vjust = 1.5, size = 5.5) +
  labs(
    x = expression("HK104 mean " * mu ~ "(" * 10^{-9} * ")"),
    y = expression("PB800 mean " * mu ~ "(" * 10^{-9} * ")")
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text  = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13)
  )

save_fig(fig_corr, "Fig_correlation_HK104_PB800", 6.5, 5.5)

# -- Fig 3: Top 15 most mutable motifs per strain --
top_motifs <- mut_counts %>%
  filter(!is.na(mean_rate), mean_rate > 0) %>%
  group_by(Strain) %>%
  slice_max(mean_rate, n = 15) %>%
  ungroup()

fig_top <- top_motifs %>%
  ggplot(aes(x = reorder(label, mean_rate),
             y = mean_rate * 1e9, fill = sub_type)) +
  geom_col(width = 0.75, color = NA) +
  geom_errorbar(
    aes(ymin = pmax(0, (mean_rate - sem_rate) * 1e9),
        ymax = (mean_rate + sem_rate) * 1e9),
    width = 0.3, linewidth = 0.4, colour = "grey30"
  ) +
  scale_fill_manual(values = SUB_COLORS, name = "Substitution") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  coord_flip() +
  facet_wrap(~Strain, nrow = 1, scales = "free") +
  labs(x = "Trinucleotide context", y = ylab_rate) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "right",
    strip.background = element_rect(fill = "grey92"),
    strip.text  = element_text(size = 13, face = "bold"),
    axis.text   = element_text(size = 11, colour = "black"),
    axis.title  = element_text(size = 13)
  )

save_fig(fig_top, "Fig_top_motifs", 11, 6)

# -- Fig 4: Nucleotide-colored context bar --
NUC_COLORS <- c(A = "#00BB00", C = "#0000BB", G = "#E47E00", T = "#BB0000")

fig_nuc_color <- mut_counts %>%
  filter(!is.na(mean_rate)) %>%
  arrange(sub_type, x, z) %>%
  mutate(label    = factor(label, levels = unique(label)),
         ref_base = substr(context_pyr, 2, 2)) %>%
  ggplot(aes(x = label, y = mean_rate * 1e9, fill = ref_base)) +
  geom_col(width = 0.8, color = NA) +
  geom_errorbar(
    aes(ymin = pmax(0, (mean_rate - sem_rate) * 1e9),
        ymax = (mean_rate + sem_rate) * 1e9),
    width = 0.4, linewidth = 0.25, colour = "grey40"
  ) +
  scale_fill_manual(values = NUC_COLORS, name = "Reference\nnucleotide") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  facet_grid(Strain ~ sub_type, scales = "free_x", space = "free_x") +
  labs(x = "Trinucleotide context", y = ylab_rate) +
  theme_tri() +
  theme(legend.position = "right",
        strip.text.y = element_text(angle = 0, size = 13))

save_fig(fig_nuc_color, "Fig_trinucleotide_nuc_colors", 16, 6.5)

# =============================================================================
# STEP 7: SUMMARY
# =============================================================================
cat("\n=== KEY MOTIFS SUMMARY ===\n")
cat("Rajaei key motifs (5'-tT/Aa-3') in C. briggsae:\n")
print(mut_counts %>%
        filter(context_pyr %in% c("TTA", "TAA") |
                 grepl("t\\[T>[ACG]\\]a", label)) %>%
        select(Strain, context_pyr, label, n_mut, mean_rate, sem_rate) %>%
        mutate(mean_1e9 = mean_rate * 1e9, sem_1e9 = sem_rate * 1e9) %>%
        arrange(Strain, desc(mean_1e9)))

cat("\nTop 5 per strain (mean +/- SEM):\n")
print(mut_counts %>%
        filter(!is.na(mean_rate), mean_rate > 0) %>%
        group_by(Strain) %>%
        slice_max(mean_rate, n = 5) %>%
        select(Strain, label, sub_type, n_mut, mean_rate, sem_rate) %>%
        mutate(mean_1e9 = mean_rate * 1e9, sem_1e9 = sem_rate * 1e9))

cat(sprintf("\nCorrelation (Pearson r) = %.4f\n", r_val))
cat("  Rajaei N2 vs PB306: r > 0.7\n")
cat("\nDone.\n")
