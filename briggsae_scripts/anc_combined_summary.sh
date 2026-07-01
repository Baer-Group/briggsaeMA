#!/usr/bin/env bash
#SBATCH --job-name=anc_combined_summary
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=anc_combined_summary_%j.out
#SBATCH --account=baer --qos=baer

# =============================================================================
# Combined SR + LR Ancestor Variant Count Summary
# =============================================================================
# CHANGES FROM PREVIOUS VERSION:
#   1. Indel sites with |len(ALT) - len(REF)| > 20 are excluded
#   2. Single output CSV with one column per tag instead of 4 separate files
#
# OUTPUT:
#   /blue/baer/m.rifat/briggsae_ancestors/combined_summary/
#     combined_count_summary.csv   — CATEGORY, snp.DP3, snp.DP10, indel.DP3, indel.DP10
# =============================================================================

set -euo pipefail

LR_DIR="/blue/baer/m.rifat/briggsae_ancestors/anc_lr_out"
SR_DIR="/blue/baer/m.rifat/briggsae_ancestors/anc_sr_out"
OUT_DIR="/blue/baer/m.rifat/briggsae_ancestors/combined_summary"

mkdir -p "$OUT_DIR"

module load python/3.10

python3 << 'PYEOF'
import os, csv
from collections import defaultdict

lr_dir  = "/blue/baer/m.rifat/briggsae_ancestors/anc_lr_out"
sr_dir  = "/blue/baer/m.rifat/briggsae_ancestors/anc_sr_out"
out_dir = "/blue/baer/m.rifat/briggsae_ancestors/combined_summary"

# ==========================================================================
# Tag pairs: (lr_filename_tag, sr_filename_tag, output_tag)
# ==========================================================================
tag_pairs = [
    ("snps.DP3",    "snp.3x",    "snp.DP3"),
    ("snps.DP10",   "snp.10x",   "snp.DP10"),
    ("indels.DP3",  "indel.3x",  "indel.DP3"),
    ("indels.DP10", "indel.10x", "indel.DP10"),
]

out_tags = [t[2] for t in tag_pairs]

# ==========================================================================
# Helper: parse sample name → (ancestor_prefix, sample_type, replicate)
# ==========================================================================
def parse_sample(name):
    first = name.split('_')[0]
    if first.upper().startswith('HK'):
        anc = 'HK'
    elif first.upper().startswith('PB'):
        anc = 'PB'
    else:
        anc = first

    if 'hifi_reads' in name or 'hifi' in name.lower() or 'G0' in name:
        return anc, 'LR', 'LR'
    else:
        base = name.split('.bp')[0] if '.bp' in name else name
        rep_token = base.split('_')[-1]
        rep = 'R1' if rep_token.endswith('1') else 'R2'
        return anc, 'SR', rep

def make_label(anc, stype, rep):
    return f"{anc}_LR" if stype == 'LR' else f"{anc}_SR_{rep}"

# ==========================================================================
# Build all 63 exact-pattern row definitions
# ==========================================================================
def build_row_defs(all_labels):
    hk_lr = next(l for l in all_labels if l == 'HK_LR')
    pb_lr = next(l for l in all_labels if l == 'PB_LR')
    hk_r1 = next(l for l in all_labels if l == 'HK_SR_R1')
    hk_r2 = next(l for l in all_labels if l == 'HK_SR_R2')
    pb_r1 = next(l for l in all_labels if l == 'PB_SR_R1')
    pb_r2 = next(l for l in all_labels if l == 'PB_SR_R2')
    S = [hk_lr, hk_r1, hk_r2, pb_lr, pb_r1, pb_r2]

    rows = []
    def add_section(title): rows.append((title, None))
    def add_row(label, pattern): rows.append((label, frozenset(pattern)))

    add_section("=== SIZE 6: All 6 samples ===")
    add_row("All 6 (HK_LR+SR_R1+R2, PB_LR+SR_R1+R2)", S)

    add_section("=== SIZE 5: 5 of 6 samples ===")
    add_section("--- HK complete, PB missing one ---")
    add_row("HK all 3 + PB_LR + PB_SR_R1  (PB_SR_R2 absent)", [hk_lr,hk_r1,hk_r2,pb_lr,pb_r1])
    add_row("HK all 3 + PB_LR + PB_SR_R2  (PB_SR_R1 absent)", [hk_lr,hk_r1,hk_r2,pb_lr,pb_r2])
    add_row("HK all 3 + PB_SR_R1 + PB_SR_R2  (PB_LR absent)", [hk_lr,hk_r1,hk_r2,pb_r1,pb_r2])
    add_section("--- PB complete, HK missing one ---")
    add_row("PB all 3 + HK_LR + HK_SR_R1  (HK_SR_R2 absent)", [pb_lr,pb_r1,pb_r2,hk_lr,hk_r1])
    add_row("PB all 3 + HK_LR + HK_SR_R2  (HK_SR_R1 absent)", [pb_lr,pb_r1,pb_r2,hk_lr,hk_r2])
    add_row("PB all 3 + HK_SR_R1 + HK_SR_R2  (HK_LR absent)", [pb_lr,pb_r1,pb_r2,hk_r1,hk_r2])

    add_section("=== SIZE 4: 4 of 6 samples ===")
    add_section("--- Both ancestors fully represented (2 from each) ---")
    add_row("HK_LR+SR_R1 + PB_LR+SR_R1  (HK_R2, PB_R2 absent)", [hk_lr,hk_r1,pb_lr,pb_r1])
    add_row("HK_LR+SR_R1 + PB_LR+SR_R2  (HK_R2, PB_R1 absent)", [hk_lr,hk_r1,pb_lr,pb_r2])
    add_row("HK_LR+SR_R1 + PB_SR_R1+R2  (HK_R2, PB_LR absent)", [hk_lr,hk_r1,pb_r1,pb_r2])
    add_row("HK_LR+SR_R2 + PB_LR+SR_R1  (HK_R1, PB_R2 absent)", [hk_lr,hk_r2,pb_lr,pb_r1])
    add_row("HK_LR+SR_R2 + PB_LR+SR_R2  (HK_R1, PB_R1 absent)", [hk_lr,hk_r2,pb_lr,pb_r2])
    add_row("HK_LR+SR_R2 + PB_SR_R1+R2  (HK_R1, PB_LR absent)", [hk_lr,hk_r2,pb_r1,pb_r2])
    add_row("HK_SR_R1+R2 + PB_LR+SR_R1  (HK_LR, PB_R2 absent)", [hk_r1,hk_r2,pb_lr,pb_r1])
    add_row("HK_SR_R1+R2 + PB_LR+SR_R2  (HK_LR, PB_R1 absent)", [hk_r1,hk_r2,pb_lr,pb_r2])
    add_row("HK_SR_R1+R2 + PB_SR_R1+R2  (both LR absent)",       [hk_r1,hk_r2,pb_r1,pb_r2])
    add_section("--- All 3 from one ancestor + 1 from other ---")
    add_row("HK all 3 + PB_LR only  (PB_SR absent)",        [hk_lr,hk_r1,hk_r2,pb_lr])
    add_row("HK all 3 + PB_SR_R1 only  (PB_LR+R2 absent)",  [hk_lr,hk_r1,hk_r2,pb_r1])
    add_row("HK all 3 + PB_SR_R2 only  (PB_LR+R1 absent)",  [hk_lr,hk_r1,hk_r2,pb_r2])
    add_row("PB all 3 + HK_LR only  (HK_SR absent)",         [pb_lr,pb_r1,pb_r2,hk_lr])
    add_row("PB all 3 + HK_SR_R1 only  (HK_LR+R2 absent)",  [pb_lr,pb_r1,pb_r2,hk_r1])
    add_row("PB all 3 + HK_SR_R2 only  (HK_LR+R1 absent)",  [pb_lr,pb_r1,pb_r2,hk_r2])

    add_section("=== SIZE 3: 3 of 6 samples ===")
    add_section("--- All 3 from one ancestor only ---")
    add_row("HK all 3 only  (PB fully absent)", [hk_lr,hk_r1,hk_r2])
    add_row("PB all 3 only  (HK fully absent)", [pb_lr,pb_r1,pb_r2])
    add_section("--- LR + one SR of same ancestor only ---")
    add_row("HK_LR + HK_SR_R1 only", [hk_lr,hk_r1])
    add_row("HK_LR + HK_SR_R2 only", [hk_lr,hk_r2])
    add_row("PB_LR + PB_SR_R1 only", [pb_lr,pb_r1])
    add_row("PB_LR + PB_SR_R2 only", [pb_lr,pb_r2])
    add_section("--- Both SR of same ancestor only (no LR) ---")
    add_row("HK_SR_R1 + HK_SR_R2 only  (HK_LR absent)", [hk_r1,hk_r2])
    add_row("PB_SR_R1 + PB_SR_R2 only  (PB_LR absent)", [pb_r1,pb_r2])
    add_section("--- Both LR + one SR (cross or same ancestor) ---")
    add_row("HK_LR + PB_LR + HK_SR_R1 only", [hk_lr,pb_lr,hk_r1])
    add_row("HK_LR + PB_LR + HK_SR_R2 only", [hk_lr,pb_lr,hk_r2])
    add_row("HK_LR + PB_LR + PB_SR_R1 only", [hk_lr,pb_lr,pb_r1])
    add_row("HK_LR + PB_LR + PB_SR_R2 only", [hk_lr,pb_lr,pb_r2])
    add_section("--- One LR + two SR from different ancestors ---")
    add_row("HK_LR + HK_SR_R1 + PB_SR_R1 only", [hk_lr,hk_r1,pb_r1])
    add_row("HK_LR + HK_SR_R1 + PB_SR_R2 only", [hk_lr,hk_r1,pb_r2])
    add_row("HK_LR + HK_SR_R2 + PB_SR_R1 only", [hk_lr,hk_r2,pb_r1])
    add_row("HK_LR + HK_SR_R2 + PB_SR_R2 only", [hk_lr,hk_r2,pb_r2])
    add_row("PB_LR + PB_SR_R1 + HK_SR_R1 only", [pb_lr,pb_r1,hk_r1])
    add_row("PB_LR + PB_SR_R1 + HK_SR_R2 only", [pb_lr,pb_r1,hk_r2])
    add_row("PB_LR + PB_SR_R2 + HK_SR_R1 only", [pb_lr,pb_r2,hk_r1])
    add_row("PB_LR + PB_SR_R2 + HK_SR_R2 only", [pb_lr,pb_r2,hk_r2])
    add_section("--- Cross-ancestor SR only (3 SR, no LR) ---")
    add_row("HK_SR_R1 + HK_SR_R2 + PB_SR_R1 only", [hk_r1,hk_r2,pb_r1])
    add_row("HK_SR_R1 + HK_SR_R2 + PB_SR_R2 only", [hk_r1,hk_r2,pb_r2])
    add_row("HK_SR_R1 + PB_SR_R1 + PB_SR_R2 only", [hk_r1,pb_r1,pb_r2])
    add_row("HK_SR_R2 + PB_SR_R1 + PB_SR_R2 only", [hk_r2,pb_r1,pb_r2])

    add_section("=== SIZE 2: 2 of 6 samples ===")
    add_section("--- Both LR only ---")
    add_row("HK_LR + PB_LR only  (no SR)", [hk_lr,pb_lr])
    add_section("--- LR + SR cross ancestor ---")
    add_row("HK_LR + PB_SR_R1 only", [hk_lr,pb_r1])
    add_row("HK_LR + PB_SR_R2 only", [hk_lr,pb_r2])
    add_row("PB_LR + HK_SR_R1 only", [pb_lr,hk_r1])
    add_row("PB_LR + HK_SR_R2 only", [pb_lr,hk_r2])
    add_section("--- SR cross ancestor (one each, no LR) ---")
    add_row("HK_SR_R1 + PB_SR_R1 only", [hk_r1,pb_r1])
    add_row("HK_SR_R1 + PB_SR_R2 only", [hk_r1,pb_r2])
    add_row("HK_SR_R2 + PB_SR_R1 only", [hk_r2,pb_r1])
    add_row("HK_SR_R2 + PB_SR_R2 only", [hk_r2,pb_r2])

    add_section("=== SIZE 1: Single sample only ===")
    add_row("HK_LR only",    [hk_lr])
    add_row("HK_SR_R1 only", [hk_r1])
    add_row("HK_SR_R2 only", [hk_r2])
    add_row("PB_LR only",    [pb_lr])
    add_row("PB_SR_R1 only", [pb_r1])
    add_row("PB_SR_R2 only", [pb_r2])

    add_section("=== GRAND TOTAL per sample (presence irrespective of others) ===")
    rows.append(("HK_LR total",    'grand', hk_lr))
    rows.append(("HK_SR_R1 total", 'grand', hk_r1))
    rows.append(("HK_SR_R2 total", 'grand', hk_r2))
    rows.append(("PB_LR total",    'grand', pb_lr))
    rows.append(("PB_SR_R1 total", 'grand', pb_r1))
    rows.append(("PB_SR_R2 total", 'grand', pb_r2))

    return rows

# ==========================================================================
# Load CSV keyed by (CHROM, POS, REF, ALT)
# For indel tags: skip sites where |len(ALT) - len(REF)| > 20
# ==========================================================================
def load_csv_keyed(fpath, is_indel=False):
    result = {}
    skipped = 0
    with open(fpath, 'r') as f:
        reader = csv.DictReader(f)
        raw_samples = [c for c in reader.fieldnames
                       if c not in ('CHROM','POS','REF','ALT')]
        labels = [make_label(*parse_sample(raw)) for raw in raw_samples]
        for row in reader:
            ref, alt = row['REF'], row['ALT']

            # Filter indels > 20 bp
            if is_indel and abs(len(alt) - len(ref)) > 20:
                skipped += 1
                continue

            key = (row['CHROM'], row['POS'], ref, alt)
            hom = frozenset(
                labels[i] for i, raw in enumerate(raw_samples)
                if set(row[raw].replace('|','/').split('/')) == {'1'}
            )
            if hom:
                result[key] = result.get(key, frozenset()) | hom

    if is_indel and skipped > 0:
        print(f"    Skipped {skipped} indel sites with |ALT-REF| > 20 bp")
    return result

# ==========================================================================
# Main: collect counts for all 4 tags, then write one combined CSV
# ==========================================================================
all_counts   = {}   # out_tag -> site_counts dict
all_row_defs = None

for lr_tag, sr_tag, out_tag in tag_pairs:
    print(f"\n{'='*60}")
    print(f"Processing: {out_tag}")
    print(f"{'='*60}")

    lr_file = os.path.join(lr_dir, f"{lr_tag}.ancestor_variants.csv")
    sr_file = os.path.join(sr_dir, f"briggsae_ancestor.bp.GVCFs.{sr_tag}.csv")

    if not os.path.exists(lr_file):
        print(f"  ERROR: LR file not found: {lr_file}"); continue
    if not os.path.exists(sr_file):
        print(f"  ERROR: SR file not found: {sr_file}"); continue

    is_indel = 'indel' in out_tag

    lr_keyed = load_csv_keyed(lr_file, is_indel=is_indel)
    sr_keyed = load_csv_keyed(sr_file, is_indel=is_indel)

    all_keys = set(lr_keyed.keys()) | set(sr_keyed.keys())

    site_counts = defaultdict(int)
    for key in all_keys:
        combined = lr_keyed.get(key, frozenset()) | sr_keyed.get(key, frozenset())
        if combined:
            site_counts[combined] += 1

    all_labels_set = set()
    for s in site_counts: all_labels_set |= s
    all_labels = sorted(all_labels_set)

    print(f"  All labels : {all_labels}")
    print(f"  Total sites after filtering: {sum(site_counts.values())}")

    all_counts[out_tag] = site_counts

    if all_row_defs is None:
        all_row_defs = build_row_defs(all_labels)

# ==========================================================================
# Write single combined CSV — columns: CATEGORY, snp.DP3, snp.DP10,
#                                       indel.DP3, indel.DP10
# ==========================================================================
out_csv = os.path.join(out_dir, "combined_count_summary.csv")

def exact_count(site_counts, pattern):
    return site_counts.get(frozenset(pattern), 0)

def grand_count(site_counts, label):
    return sum(v for k, v in site_counts.items() if label in k)

with open(out_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['CATEGORY'] + out_tags)
    for row in all_row_defs:
        if row[1] is None:
            # Section header — empty count columns
            writer.writerow([row[0]] + ['']*len(out_tags))
        elif row[1] == 'grand':
            writer.writerow([row[0]] + [grand_count(all_counts[t], row[2]) for t in out_tags])
        else:
            writer.writerow([row[0]] + [exact_count(all_counts[t], row[1]) for t in out_tags])

print(f"\n\n[DONE] Written: {out_csv}", flush=True)
PYEOF

echo ""
echo "================================================================"
echo "[DONE] Combined summary complete."
echo "Output: $OUT_DIR/combined_count_summary.csv"
echo "================================================================"
