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
# If briggsae_genome_trinucleotide_counts.csv already stores callable motif counts,
  # do not multiply by a global callable fraction again.
  mutate(callable_count = genome_count)

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
  mutate(rate = n_mut / (callable_count * GENERATIONS*CALLABLE_FRAC))

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
      legend.title     = element_text(size = 20),
      legend.text      = element_text(size = 20),
      legend.position  = "none",
      panel.spacing    = unit(0.1, "lines"),
      plot.title       = element_blank(),
      plot.subtitle    = element_text(size = 22, face = "bold")
    )
}

save_fig <- function(p, name, w, h) {
  ggsave(paste0(name, ".pdf"), p, width = w, height = h)
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s\n", name))
}

# Y-axis label: pure plotmath, no Unicode
ylab_rate <- expression("Mean " * mu ~ "(" * 10^{-9} ~ "per site per gen)")

# -- Fig 1: substitution-specific trinucleotide context rates with SEM error bars --
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

# -- Additional Rajaei/Saxena-style Fig: 64 three-bp motif rates --
# This collapses over substitution outcome and estimates one base-substitution
# mutation rate for each exact 5'-xYz-3' motif, with the mutable base Y in the
# middle position. This is the 64-motif representation used by Rajaei et al.
# (2021) and Saxena et al. (2019), rather than the 96 substitution-context
# representation used above.
cat("Computing 64 exact three-bp motif rates...\n")

tri_genome_64 <- read.csv("briggsae_genome_trinucleotide_counts.csv",
                          stringsAsFactors = FALSE) %>%
  mutate(trinucleotide = toupper(trinucleotide)) %>%
  filter(nchar(trinucleotide) == 3, grepl("^[ACGT]{3}$", trinucleotide)) %>%
  group_by(context_64 = trinucleotide) %>%
  summarise(genome_count_64 = sum(count), .groups = "drop") %>%
  # If briggsae_genome_trinucleotide_counts.csv already stores callable motif counts,
  # do not multiply by a global callable fraction again.
  mutate(callable_count_64 = genome_count_64)

motif_levels <- expand.grid(
  x = c("A", "C", "G", "T"),
  Y = c("A", "C", "G", "T"),
  z = c("A", "C", "G", "T"),
  stringsAsFactors = FALSE
) %>%
  arrange(Y, x, z) %>%
  mutate(context_64 = paste0(x, Y, z),
         label_64   = paste0(tolower(x), Y, tolower(z)))

full_grid_64 <- line_strain %>%
  cross_join(tri_genome_64 %>% select(context_64, callable_count_64)) %>%
  left_join(
    muts %>% count(Sample, context_64 = context, name = "n_mut_64"),
    by = c("Sample", "context_64")
  ) %>%
  replace_na(list(n_mut_64 = 0)) %>%
  mutate(rate_64 = n_mut_64 / (callable_count_64 * GENERATIONS*CALLABLE_FRAC))

mut_counts_64 <- full_grid_64 %>%
  group_by(Strain, context_64) %>%
  summarise(
    n_lines      = n(),
    n_mut        = sum(n_mut_64),
    mean_rate    = mean(rate_64),
    sd_rate      = sd(rate_64),
    sem_rate     = sd(rate_64) / sqrt(n()),
    .groups      = "drop"
  ) %>%
  left_join(tri_genome_64, by = "context_64") %>%
  left_join(motif_levels, by = "context_64") %>%
  mutate(
    focal_base = substr(context_64, 2, 2),
    label_64   = factor(label_64, levels = motif_levels$label_64)
  )

write.csv(mut_counts_64, "trinucleotide_64motif_rates.csv", row.names = FALSE)
cat("Saved: trinucleotide_64motif_rates.csv\n")

NUC_COLORS <- c(A = "#00BB00", C = "#0000BB", G = "#E47E00", T = "#BB0000")

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
      legend.title     = element_text(size = 26),
      legend.text      = element_text(size = 26),
      panel.spacing    = unit(0.12, "lines"),
      plot.title       = element_blank(),
      plot.subtitle    = element_text(size = 28, face = "bold")
    )
}

plot_64motif <- function(strain_name) {
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
    scale_fill_manual(values = NUC_COLORS, name = "Middle base") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_grid(. ~ focal_base, scales = "free_x", space = "free_x") +
    labs(x = "64 three-bp motif (5'-xYz-3'; Y = middle base)",
         y = ylab_rate,
         subtitle = strain_name) +
    theme_tri64()
}

p64_hk <- plot_64motif("HK104")
p64_pb <- plot_64motif("PB800")
fig_64_combined <- p64_hk / p64_pb +
  plot_layout(guides = "collect") &
  theme(legend.position = "top",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))

save_fig(fig_64_combined, "Fig_trinucleotide_64motif_combined", 16, 9.5)

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
  theme_bw(base_size = 18) +
  theme(
    legend.position = "right",
    axis.text  = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
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
  theme_bw(base_size = 18) +
  theme(
    legend.position  = "right",
    strip.background = element_rect(fill = "grey92"),
    strip.text  = element_text(size = 18, face = "bold"),
    axis.text   = element_text(size = 18, colour = "black"),
    axis.title  = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
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
        strip.text.y = element_text(angle = 0, size = 18))

save_fig(fig_nuc_color, "Fig_trinucleotide_nuc_colors", 16, 6.5)


# =============================================================================
# STEP 6b: STRAIN-DIFFERENCE BAR CHARTS  (HK104 - PB800)
# -----------------------------------------------------------------------------
# Diverging bar chart: each bar = (HK104 rate - PB800 rate) per context.
#   - Bars ABOVE zero  -> HK104 higher  (blue,  HK_DIFF_BLUE)
#   - Bars BELOW zero  -> PB800 higher  (orange, PB_DIFF_ORANGE)
# Context ORDER IS PRESERVED (NOT re-sorted by magnitude): the x-axis order
# and faceting match the existing spectrum figures exactly, so this difference
# plot lines up positionally with Fig_trinucleotide_combined (96-channel) and
# Fig_trinucleotide_64motif_combined (64-motif).
# =============================================================================
cat("\n--- Step 6b: strain-difference bar charts (HK104 - PB800) ---\n")

HK_DIFF_BLUE   <- "#1f4e8f"   # HK104 higher
PB_DIFF_ORANGE <- "#e67e22"   # PB800 higher

# How many of the most-divergent contexts to label with their value.
# Set to 0 to disable all value labels (cleanest); 8-12 highlights the extremes.
N_LABEL_EXTREMES <- 10

# Shared y-axis label (plotmath, no Unicode)
ylab_diff <- expression(Delta * mu~"(HK104 - PB800, in"~10^{-9}*")")

# -----------------------------------------------------------------------------
# (A) 96-channel difference  (substitution-specific context)
#     Order preserved from plot_trinucleotide(): arrange(sub_type, x, z)
# -----------------------------------------------------------------------------
order_96 <- mut_counts %>%
  distinct(sub_type, x, z, label) %>%
  arrange(sub_type, x, z) %>%
  pull(label)

diff_96 <- mut_counts %>%
  select(Strain, label, sub_type, x, z, mean_rate) %>%
  pivot_wider(names_from = Strain, values_from = mean_rate, values_fill = 0) %>%
  mutate(
    diff   = (HK104 - PB800) * 1e9,            # Delta-mu in 1e-9 units
    higher = ifelse(diff >= 0, "HK104", "PB800"),
    label  = factor(label, levels = order_96)
  ) %>%
  arrange(label) %>%
  mutate(abs_diff = abs(diff),
         lab_val  = ifelse(rank(-abs_diff, ties.method = "first") <= N_LABEL_EXTREMES,
                           sprintf("%+.2f", diff), NA_character_))

ymax_96 <- max(abs(diff_96$diff), na.rm = TRUE) * 1.20

fig_diff_96 <- ggplot(diff_96, aes(x = label, y = diff, fill = higher)) +
  geom_col(width = 0.8, color = NA) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c(HK104 = HK_DIFF_BLUE, PB800 = PB_DIFF_ORANGE),
                    name = "Strain") +
  scale_y_continuous(limits = c(-ymax_96, ymax_96)) +
  facet_grid(. ~ sub_type, scales = "free_x", space = "free_x") +
  labs(x = "Base Substitution-Local contexts",
       y = ylab_diff) +
  theme_tri() +
  theme(legend.position = "top",
        legend.title = element_text(size = 18),
        legend.text  = element_text(size = 18))

save_fig(fig_diff_96, "Fig_strain_diff_96channel", 16, 8)

# -----------------------------------------------------------------------------
# (B) 64-motif difference  (exact 5'-xYz-3' motif)
#     Order preserved from plot_64motif(): motif_levels$label_64 factor order
# -----------------------------------------------------------------------------
diff_64 <- mut_counts_64 %>%
  select(Strain, context_64, label_64, focal_base, mean_rate) %>%
  pivot_wider(names_from = Strain, values_from = mean_rate, values_fill = 0) %>%
  mutate(
    diff     = (HK104 - PB800) * 1e9,
    higher   = ifelse(diff >= 0, "HK104", "PB800"),
    label_64 = factor(label_64, levels = motif_levels$label_64)
  ) %>%
  arrange(label_64) %>%
  mutate(abs_diff = abs(diff),
         lab_val  = ifelse(rank(-abs_diff, ties.method = "first") <= N_LABEL_EXTREMES,
                           sprintf("%+.2f", diff), NA_character_))

ymax_64 <- max(abs(diff_64$diff), na.rm = TRUE) * 1.20

fig_diff_64 <- ggplot(diff_64, aes(x = label_64, y = diff, fill = higher)) +
  geom_col(width = 0.8, color = NA) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c(HK104 = HK_DIFF_BLUE, PB800 = PB_DIFF_ORANGE),
                    name = "Higher strain") +
  scale_y_continuous(limits = c(-ymax_64, ymax_64)) +
  facet_grid(. ~ focal_base, scales = "free_x", space = "free_x") +
  labs(x = "64 three-bp motif (5'-xYz-3'; canonical order - position preserved)",
       y = ylab_diff) +
  theme_tri64() +
  theme(legend.position = "top",
        legend.title = element_text(size = 18),
        legend.text  = element_text(size = 18))

save_fig(fig_diff_64, "Fig_strain_diff_64motif", 16, 8)

# -----------------------------------------------------------------------------
# Console summary: most strain-divergent contexts
# -----------------------------------------------------------------------------
cat("\nTop strain-divergent 96-channel contexts (|HK104 - PB800|):\n")
print(diff_96 %>% arrange(desc(abs_diff)) %>%
        transmute(label, sub_type,
                  HK104_1e9 = HK104 * 1e9, PB800_1e9 = PB800 * 1e9, diff_1e9 = diff) %>%
        head(N_LABEL_EXTREMES))

cat("\nTop strain-divergent 64-motif contexts (|HK104 - PB800|):\n")
print(diff_64 %>% arrange(desc(abs_diff)) %>%
        transmute(label_64, focal_base,
                  HK104_1e9 = HK104 * 1e9, PB800_1e9 = PB800 * 1e9, diff_1e9 = diff) %>%
        head(N_LABEL_EXTREMES))

cat("\nSaved: Fig_strain_diff_96channel.pdf/.png, Fig_strain_diff_64motif.pdf/.png\n")

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

