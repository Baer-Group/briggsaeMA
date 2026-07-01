#!/usr/bin/env Rscript
# sv_plots.R — v2 (cleaned-up plots, no titles)
# -----------------------------------------------------------------------------
# Generates the two SV figures from de_novo_hom_carriers.tsv:
#   1) fig_sv_per_line_hom.png — stacked bar of per-line counts by SVTYPE,
#                                with strain-group labels and a divider between
#                                HK104 and PB800 panels (no plot title).
#   2) fig_sv_size_hom.png     — faceted size distribution by SVTYPE × strain
#                                (no plot title).
#
# Usage:
#   Rscript sv_plots.R [input.tsv] [out_dir]
#   # or, in interactive R:
#   source("sv_plots.R")
#
# Defaults:
#   input.tsv = de_novo_hom_carriers.tsv in the working directory
#   out_dir   = working directory
#
# Input TSV columns (header expected):
#   chrom  pos  svtype  svlen  id  sample  gt  dp
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

# ---- Configuration ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_tsv <- if (length(args) >= 1) args[1] else "de_novo_hom_carriers.tsv"
out_dir   <- if (length(args) >= 2) args[2] else "."

svtype_colors <- c(
  DEL = "#1f77b4",
  INS = "#ff7f0e",
  DUP = "#2ca02c",
  INV = "#d62728",
  BND = "#9467bd"
)
strain_colors <- c(
  HK104 = "#162A6F",   # blue4-like
  PB800 = "#FF7F50"    # coral
)

# ---- Read -------------------------------------------------------------------
if (!file.exists(input_tsv)) {
  stop("Input file not found: ", input_tsv,
       "\nPass the path as the first argument: Rscript sv_plots.R /path/to/file.tsv")
}

dn <- read_tsv(input_tsv, show_col_types = FALSE)

required <- c("chrom","pos","svtype","svlen","sample","gt","dp")
missing  <- setdiff(required, names(dn))
if (length(missing) > 0) {
  stop("Missing expected column(s): ", paste(missing, collapse = ", "))
}

# Derive strain, clean line name, parse SV length
dn <- dn |>
  mutate(
    line       = sub("\\.hifi_reads$", "", sample),
    line_clean = sub("_", "", line),                # HK_205 -> HK205 for display
    strain     = if_else(startsWith(line, "HK"), "HK104", "PB800"),
    svtype     = factor(svtype, levels = c("DEL","INS","DUP","INV","BND")),
    svlen_num  = suppressWarnings(as.numeric(svlen)),
    abs_svlen  = abs(svlen_num)
  )

cat("Total calls:", nrow(dn), "\n")
cat("\nBy strain:\n");  print(table(dn$strain))
cat("\nBy SVTYPE:\n");  print(table(dn$svtype, useNA = "ifany"))

# ---- Per-line × SVTYPE counts -----------------------------------------------
counts_wide <- dn |>
  count(strain, line, line_clean, svtype, .drop = FALSE) |>
  pivot_wider(names_from = svtype, values_from = n, values_fill = 0) |>
  mutate(total = DEL + INS + DUP + INV + BND) |>
  arrange(strain, line)

cat("\n=== Per-line counts ===\n")
print(counts_wide)

strain_means <- counts_wide |>
  group_by(strain) |>
  summarise(mean_per_line = mean(total),
            n_lines       = n(),
            total         = sum(total),
            .groups       = "drop")
cat("\n=== Per-strain summary ===\n")
print(strain_means)
ratio <- strain_means$mean_per_line[strain_means$strain == "HK104"] /
         strain_means$mean_per_line[strain_means$strain == "PB800"]
cat(sprintf("\nHK104 / PB800 per-line ratio: %.3fx\n", ratio))

# Stable factor ordering across HK then PB
ordered_lines <- counts_wide$line_clean

counts_long <- dn |>
  count(strain, line_clean, svtype, .drop = FALSE) |>
  mutate(line_clean = factor(line_clean, levels = ordered_lines))

totals_per_line <- counts_wide |>
  mutate(line_clean = factor(line_clean, levels = ordered_lines))

# Geometry for divider and strain-group labels
n_hk      <- sum(counts_wide$strain == "HK104")
n_pb      <- sum(counts_wide$strain == "PB800")
divider_x <- n_hk + 0.5
hk_center <- (1 + n_hk) / 2
pb_center <- n_hk + (1 + n_pb) / 2

# Group labels sit just below the x-axis tick labels.
# coord_cartesian(clip="off") allows drawing outside the panel;
# plot.margin gives room.
y_max       <- max(counts_wide$total)
y_label_pos <- -y_max * 0.10

# ---- Figure 1: per-line stacked bar -----------------------------------------
p1 <- ggplot(counts_long, aes(x = line_clean, y = n, fill = svtype)) +
  geom_col(color = "white", linewidth = 0.4) +
  geom_text(data = totals_per_line,
            aes(x = line_clean, y = total, label = total),
            inherit.aes = FALSE,
            vjust = -0.4, fontface = "bold", size = 4.5) +
  geom_vline(xintercept = divider_x,
             linetype = "dashed", color = "grey50",
             linewidth = 0.5, alpha = 0.6) +
  annotate("text",
           x = hk_center, y = y_label_pos, label = "HK104",
           color = strain_colors[["HK104"]], fontface = "bold", size = 5.5,
           hjust = 0.5, vjust = 1) +
  annotate("text",
           x = pb_center, y = y_label_pos, label = "PB800",
           color = strain_colors[["PB800"]], fontface = "bold", size = 5.5,
           hjust = 0.5, vjust = 1) +
  scale_fill_manual(values = svtype_colors, name = "SVTYPE", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Singleton SV count") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12),
    axis.text.y        = element_text(size = 11),
    axis.title.y       = element_text(size = 14),
    legend.position    = "right",
    legend.title       = element_text(size = 13, face = "bold"),
    legend.text        = element_text(size = 12),
    plot.margin        = margin(t = 10, r = 12, b = 36, l = 10)
  )

out1 <- file.path(out_dir, "fig_sv_per_line_hom.png")
ggsave(out1, p1, width = 12, height = 5.5, dpi = 200)
cat("Saved:", out1, "\n")

# ---- Figure 2: size distribution by SVTYPE (no title, larger fonts) ---------
sized <- dn |>
  filter(svtype %in% c("DEL","INS","DUP","INV"), !is.na(abs_svlen))

p2 <- ggplot(sized, aes(x = abs_svlen, fill = strain)) +
  geom_histogram(bins = 25, color = "white", linewidth = 0.3,
                 position = "identity", alpha = 0.65) +
  facet_wrap(~ svtype, scales = "free_y") +
  scale_x_log10(labels = label_comma()) +
  scale_fill_manual(values = strain_colors, name = "Strain") +
  labs(x = expression("|SV length| (bp, log"[10]*" scale)"),
       y = "Count") +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 13),
    axis.text        = element_text(size = 11),
    axis.title       = element_text(size = 13),
    legend.position  = "top",
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 12)
  )

out2 <- file.path(out_dir, "fig_sv_size_hom.png")
ggsave(out2, p2, width = 11, height = 7, dpi = 200)
cat("Saved:", out2, "\n")

cat("\nDone.\n")
