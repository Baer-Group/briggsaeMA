#!/usr/bin/env bash
#SBATCH --job-name=anc_sr_variants
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --time=96:00:00
#SBATCH --output=anc_sr_variants_%j.out
#SBATCH --account=baer --qos=baer

# =============================================================================
# C. briggsae Short-Read Ancestor Variant Listing Pipeline
# =============================================================================
# OVERVIEW:
#   Joint VCF (briggsae_ancestor.bp.GVCFs.vcf) assumed to already exist.
#   For each variant type (SNP / INDEL) and coverage threshold (3x / 10x):
#       1. Split variant type (GATK SelectVariants)
#       2. Mask low-DP genotypes (GATK VariantFiltration + SelectVariants)
#       3. Keep only biallelic sites (bcftools)
#       4. AF filter: for 1/1 samples only, mask to ./. if minority AF > 10%
#          Heterozygotes kept as-is.
#       5. Keep sites with at least one 1/1 remaining
#       6. Export CSV: CHROM, POS, REF, ALT, [sample_name columns]
#       7. Generate 19-row count summary across all 4 samples:
#            Samples parsed from filename — first token = ancestor (HK104/PB800),
#            last token before .bp = replicate (L001/L002)
#
# COUNT STRUCTURE (19 rows per tag):
#   Category 1 — All 4 samples (1 row)
#   Category 2 — Both replicates of one ancestor only (2 rows)
#   Category 3 — Unique to one replicate (4 rows)
#   Category 4 — Cross-contamination: one replicate of each ancestor (4 rows)
#   Category 5 — Both of one ancestor + one replicate of other (4 rows)
#   Category 6 — Grand total per sample (4 rows)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# SECTION 1: User settings
# -----------------------------------------------------------------------------

WORKDIR="/blue/baer/m.rifat/briggsae_ancestors"
REF="20250626_c_briggsae_Feb2020.genome.fa"
BASE="briggsae_ancestor.bp.GVCFs"
OUTDIR="anc_sr_out"

DP_LEVELS=(3 10)
AF_THRESHOLD=0.10
JAVA_MEM="-Xmx39g"

# -----------------------------------------------------------------------------
# SECTION 2: Setup
# -----------------------------------------------------------------------------

module purge || true
module load gatk/
module load bcftools
module load samtools
cd "$WORKDIR"
mkdir -p "$OUTDIR"

[[ -f "$REF"    ]] || { echo "ERROR: Reference not found: $REF" >&2; exit 1; }
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

JOINT="${BASE}.vcf"
[[ -f "$JOINT" ]] || { echo "ERROR: Joint VCF not found: $JOINT" >&2; exit 1; }
echo "[INFO] Joint VCF found: $JOINT"

# -----------------------------------------------------------------------------
# SECTION 3: Split SNPs and indels
# -----------------------------------------------------------------------------
echo "[STEP 1] Splitting SNPs and indels..."

SNP="${BASE}.snp.vcf"
INDEL="${BASE}.indel.vcf"

gatk --java-options "$JAVA_MEM" SelectVariants \
    -R "$REF" -V "$JOINT" \
    --select-type-to-include SNP \
    -O "$SNP"

gatk --java-options "$JAVA_MEM" SelectVariants \
    -R "$REF" -V "$JOINT" \
    --select-type-to-include INDEL \
    -O "$INDEL"

echo "[STEP 1] Done: $SNP, $INDEL"

# -----------------------------------------------------------------------------
# SECTION 4: DP filtering function
# NOTE: writes output path to temp file to avoid capturing GATK logs in $()
# -----------------------------------------------------------------------------

run_dp_filter() {
    local INPUT="$1"
    local DP="$2"
    local TAG="$3"
    local RETFILE="$4"

    local FILTERED="${BASE}.${TAG}.vcf"
    local NOCALL="${BASE}.${TAG}.nocall.vcf"
    local CLEAN="${BASE}.${TAG}.filtered.vcf"

    echo "[STEP 2] DP${DP} filter on $INPUT..."

    gatk --java-options "$JAVA_MEM" VariantFiltration \
        -R "$REF" -V "$INPUT" \
        --genotype-filter-name "DP" \
        --genotype-filter-expression "DP < ${DP}.0" \
        -O "$FILTERED"

    gatk --java-options "$JAVA_MEM" SelectVariants \
        -V "$FILTERED" \
        --set-filtered-gt-to-nocall \
        -O "$NOCALL"

    gatk --java-options "$JAVA_MEM" SelectVariants \
        -R "$REF" -V "$NOCALL" \
        --exclude-filtered true \
        --exclude-non-variants true \
        -O "$CLEAN"

    echo "$CLEAN" > "$RETFILE"
}

run_dp_filter "$SNP"   3  "snp.3x"    .snp_3x_out
run_dp_filter "$INDEL" 3  "indel.3x"  .indel_3x_out
run_dp_filter "$SNP"   10 "snp.10x"   .snp_10x_out
run_dp_filter "$INDEL" 10 "indel.10x" .indel_10x_out

SNP_3X_FILTERED=$(cat .snp_3x_out)
INDEL_3X_FILTERED=$(cat .indel_3x_out)
SNP_10X_FILTERED=$(cat .snp_10x_out)
INDEL_10X_FILTERED=$(cat .indel_10x_out)
rm -f .snp_3x_out .indel_3x_out .snp_10x_out .indel_10x_out

echo "[STEP 2] Done."

# -----------------------------------------------------------------------------
# SECTION 5: Biallelic filter
# -----------------------------------------------------------------------------
echo "[STEP 3] Applying biallelic filter..."

SNP_3X_BI="${BASE}.snp.3x.filtered.biallelic.vcf"
INDEL_3X_BI="${BASE}.indel.3x.filtered.biallelic.vcf"
SNP_10X_BI="${BASE}.snp.10x.filtered.biallelic.vcf"
INDEL_10X_BI="${BASE}.indel.10x.filtered.biallelic.vcf"

bcftools view --max-alleles 2 "$SNP_3X_FILTERED"   -o "$SNP_3X_BI"
bcftools view --max-alleles 2 "$INDEL_3X_FILTERED"  -o "$INDEL_3X_BI"
bcftools view --max-alleles 2 "$SNP_10X_FILTERED"   -o "$SNP_10X_BI"
bcftools view --max-alleles 2 "$INDEL_10X_FILTERED" -o "$INDEL_10X_BI"

echo "[STEP 3] Done."

# Unload gatk before python steps — GATK bundles its own incomplete conda
# Python which hijacks python/3.10 and causes ModuleNotFoundError: encodings
module unload gatk
module load python/3.10

# -----------------------------------------------------------------------------
# SECTION 6: AF filter + CSV export function
# -----------------------------------------------------------------------------
# For each biallelic VCF:
#   - Mask 1/1 samples to ./. if minority AF (REF/total) > 10%
#   - Keep hets as-is
#   - Drop site if no 1/1 remains
#   - Export CSV with sample GT columns

process_vcf() {
    local INPUT_VCF="$1"
    local OUTPUT_CSV="$2"
    local TAG="$3"

python3 << PYEOF
import sys

input_vcf  = "${INPUT_VCF}"
output_csv = "${OUTPUT_CSV}"
tag        = "${TAG}"
af_thresh  = float("${AF_THRESHOLD}")

MAX_MISSING = 1   # with 4 samples, allow at most 1 missing

sites_total   = 0
sites_no_hom  = 0
sites_missing = 0
sites_kept    = 0
gt_masked     = 0

samples = []

with open(input_vcf, 'r') as fin, open(output_csv, 'w') as fout:
    for line in fin:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            fields  = line.rstrip('\n').split('\t')
            samples = fields[9:]
            fout.write('CHROM,POS,REF,ALT,' + ','.join(samples) + '\n')
            continue

        sites_total += 1
        fields = line.rstrip('\n').split('\t')
        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]

        fmt = fields[8].split(':')
        try:
            gt_idx = fmt.index('GT')
        except ValueError:
            sites_total -= 1
            continue
        try:
            ad_idx = fmt.index('AD')
        except ValueError:
            ad_idx = None

        sample_data = fields[9:]
        new_samples = []

        # AF masking: 1/1 samples only
        for s in sample_data:
            sf  = s.split(':')
            gt  = sf[gt_idx].replace('|', '/') if gt_idx < len(sf) else './.'
            als = gt.split('/')

            if set(als) == {'1'} and ad_idx is not None:
                if ad_idx < len(sf) and sf[ad_idx] not in ('.', ''):
                    try:
                        ad_vals = [int(x) for x in sf[ad_idx].split(',')]
                        total   = sum(ad_vals)
                        if total > 0 and ad_vals[0] / total > af_thresh:
                            sf[gt_idx] = './.'
                            new_samples.append(':'.join(sf))
                            gt_masked += 1
                            continue
                    except ValueError:
                        pass

            new_samples.append(s)

        # Must have at least one 1/1 after masking
        has_hom_alt = any(
            set(s.split(':')[gt_idx].replace('|', '/').split('/')) == {'1'}
            for s in new_samples
            if gt_idx < len(s.split(':'))
        )
        if not has_hom_alt:
            sites_no_hom += 1
            continue

        # At most 1 missing (with 4 samples)
        missing_count = sum(
            1 for s in new_samples
            if gt_idx < len(s.split(':')) and
               '.' in s.split(':')[gt_idx].replace('|', '/').split('/')
        )
        if missing_count > MAX_MISSING:
            sites_missing += 1
            continue

        gt_vals = [
            s.split(':')[gt_idx] if gt_idx < len(s.split(':')) else './.'
            for s in new_samples
        ]
        fout.write(','.join([chrom, pos, ref, alt] + gt_vals) + '\n')
        sites_kept += 1

print(f"  [{tag}] Filter summary:")
print(f"    Sites evaluated      : {sites_total}")
print(f"    Genotypes AF-masked  : {gt_masked}")
print(f"    No hom-alt remaining : {sites_no_hom}")
print(f"    >1 missing           : {sites_missing}")
print(f"    Sites kept           : {sites_kept}")
print(f"    Output CSV           : {output_csv}", flush=True)
PYEOF
}

echo "[STEP 4] Applying AF filter and exporting CSVs..."
process_vcf "$SNP_3X_BI"    "$OUTDIR/${BASE}.snp.3x.csv"    "snp.3x"
process_vcf "$INDEL_3X_BI"  "$OUTDIR/${BASE}.indel.3x.csv"  "indel.3x"
process_vcf "$SNP_10X_BI"   "$OUTDIR/${BASE}.snp.10x.csv"   "snp.10x"
process_vcf "$INDEL_10X_BI" "$OUTDIR/${BASE}.indel.10x.csv" "indel.10x"
echo "[STEP 4] Done."

# -----------------------------------------------------------------------------
# SECTION 7: 19-row count summary
# -----------------------------------------------------------------------------
# Sample name parsing:
#   First '_'-delimited token  = ancestor ID (e.g. HK104, PB800)
#   Last token before .bp      = replicate   (e.g. L001, L002)
# The 4 samples map to: H1, H2, P1, P2
# where H = HK104 ancestor, P = PB800 ancestor, 1/2 = L001/L002
# -----------------------------------------------------------------------------

echo ""
echo "[INFO] Generating 19-row count summary..."

python3 << PYEOF
import os, csv
from collections import defaultdict
from itertools import combinations

outdir  = "${OUTDIR}"
out_csv = os.path.join(outdir, "ancestor_count_summary.csv")

files = sorted(
    f for f in os.listdir(outdir)
    if f.endswith('.csv') and 'count' not in f
)

if not files:
    print("  No CSV files found.")
else:
    all_tags = []
    # counts[tag][frozenset_of_samples] = number of sites
    counts = {}

    for fname in files:
        tag = fname.split('.bp.GVCFs.')[1].replace('.csv', '')
        all_tags.append(tag)
        counts[tag] = defaultdict(int)

        fpath = os.path.join(outdir, fname)
        with open(fpath, 'r') as f:
            reader = csv.DictReader(f)
            samples = [c for c in reader.fieldnames
                       if c not in ('CHROM','POS','REF','ALT')]

            # Parse sample → ancestor and replicate
            # e.g. HK104_ANC_1_POOLRET91_S18_L001 → ancestor=HK104, rep=L001
            def parse_sample(s):
                parts = s.split('_')
                ancestor = parts[0]              # HK104 or PB800
                replicate = s.split('.')[-1] if '.' in s else parts[-1]
                # replicate is last underscore token before .bp
                base = s.split('.bp')[0]
                rep = base.split('_')[-1]        # L001 or L002
                return ancestor, rep

            sample_info = {s: parse_sample(s) for s in samples}

            for row in reader:
                # Which samples are 1/1 at this site
                hom_set = frozenset(
                    s for s in samples
                    if set(row[s].replace('|','/').split('/')) == {'1'}
                )
                if hom_set:
                    counts[tag][hom_set] += 1

    # Identify the two ancestors and their replicates from sample names
    # We'll use the first CSV to get sample names
    first_csv = os.path.join(outdir, files[0])
    with open(first_csv) as f:
        reader   = csv.DictReader(f)
        samples  = [c for c in reader.fieldnames
                    if c not in ('CHROM','POS','REF','ALT')]

    # Group samples by ancestor
    ancestor_groups = defaultdict(list)
    for s in samples:
        anc = s.split('_')[0]
        ancestor_groups[anc].append(s)

    ancestors = sorted(ancestor_groups.keys())
    anc1, anc2 = ancestors[0], ancestors[1]  # e.g. HK104, PB800

    a1r1, a1r2 = sorted(ancestor_groups[anc1])
    a2r1, a2r2 = sorted(ancestor_groups[anc2])

    # All 4 sample labels for display
    lbl = {a1r1: f"{anc1}_R1", a1r2: f"{anc1}_R2",
           a2r1: f"{anc2}_R1", a2r2: f"{anc2}_R2"}

    def fmt(s_set):
        return '+'.join(lbl[s] for s in sorted(s_set, key=lambda x: lbl[x]))

    def count_pattern(tag, pattern_set):
        """Count sites where exactly the given set of samples are 1/1."""
        return counts[tag].get(frozenset(pattern_set), 0)

    def count_presence(tag, required_set):
        """Count sites where at least the required samples are 1/1
           (may also be 1/1 in others)."""
        return sum(
            v for k, v in counts[tag].items()
            if frozenset(required_set).issubset(k)
        )

    def count_grand(tag, sample):
        """Count all sites where this sample is 1/1 regardless of others."""
        return sum(v for k, v in counts[tag].items() if sample in k)

    # Build 19-row structure
    row_defs = [
        # Category 1: All 4
        ("--- Category 1: Shared in all 4 ---",            None),
        (f"All 4 ({anc1}_R1+R2, {anc2}_R1+R2)",           lambda t: count_pattern(t, [a1r1,a1r2,a2r1,a2r2])),

        # Category 2: Both replicates of one ancestor only
        ("--- Category 2: Ancestor-specific (both reps) ---", None),
        (f"{anc1} both reps only",                         lambda t: count_pattern(t, [a1r1,a1r2])),
        (f"{anc2} both reps only",                         lambda t: count_pattern(t, [a2r1,a2r2])),

        # Category 3: Unique to one replicate
        ("--- Category 3: Unique to one replicate ---",    None),
        (f"{anc1}_R1 only",                                lambda t: count_pattern(t, [a1r1])),
        (f"{anc1}_R2 only",                                lambda t: count_pattern(t, [a1r2])),
        (f"{anc2}_R1 only",                                lambda t: count_pattern(t, [a2r1])),
        (f"{anc2}_R2 only",                                lambda t: count_pattern(t, [a2r2])),

        # Category 4: Cross-contamination (one rep each ancestor)
        ("--- Category 4: Cross (one rep each ancestor) ---", None),
        (f"{anc1}_R1 + {anc2}_R1 only",                   lambda t: count_pattern(t, [a1r1,a2r1])),
        (f"{anc1}_R1 + {anc2}_R2 only",                   lambda t: count_pattern(t, [a1r1,a2r2])),
        (f"{anc1}_R2 + {anc2}_R1 only",                   lambda t: count_pattern(t, [a1r2,a2r1])),
        (f"{anc1}_R2 + {anc2}_R2 only",                   lambda t: count_pattern(t, [a1r2,a2r2])),

        # Category 5: Both reps of one + one rep of other
        ("--- Category 5: Both reps one + one rep other ---", None),
        (f"{anc1} both + {anc2}_R1 only",                 lambda t: count_pattern(t, [a1r1,a1r2,a2r1])),
        (f"{anc1} both + {anc2}_R2 only",                 lambda t: count_pattern(t, [a1r1,a1r2,a2r2])),
        (f"{anc2} both + {anc1}_R1 only",                 lambda t: count_pattern(t, [a2r1,a2r2,a1r1])),
        (f"{anc2} both + {anc1}_R2 only",                 lambda t: count_pattern(t, [a2r1,a2r2,a1r2])),

        # Category 6: Grand total per sample
        ("--- Category 6: Grand total per sample ---",     None),
        (f"{anc1}_R1 total",                               lambda t: count_grand(t, a1r1)),
        (f"{anc1}_R2 total",                               lambda t: count_grand(t, a1r2)),
        (f"{anc2}_R1 total",                               lambda t: count_grand(t, a2r1)),
        (f"{anc2}_R2 total",                               lambda t: count_grand(t, a2r2)),
    ]

    # Print to console
    col_w  = 20
    header = f"{'CATEGORY':<55}" + "".join(f"{t:>{col_w}}" for t in all_tags)
    sep    = "-" * (55 + col_w * len(all_tags))
    print("")
    print("  Ancestor variant count summary (19 rows)")
    print("  " + sep)
    print("  " + header)
    print("  " + sep)
    for label, fn in row_defs:
        if fn is None:
            print("  " + label)
        else:
            vals = "".join(f"{fn(t):>{col_w}}" for t in all_tags)
            print(f"  {label:<55}" + vals)
    print("  " + sep)

    # Write CSV
    data_rows = [(label, fn) for label, fn in row_defs if fn is not None]
    section_rows = [(label, fn) for label, fn in row_defs]

    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['CATEGORY'] + all_tags)
        for label, fn in row_defs:
            if fn is None:
                writer.writerow([label] + ['']*len(all_tags))
            else:
                writer.writerow([label] + [fn(t) for t in all_tags])

    print(f"\n  Written: {out_csv}", flush=True)
PYEOF

# -----------------------------------------------------------------------------
# SECTION 8: Cleanup intermediates
# -----------------------------------------------------------------------------
echo ""
echo "[CLEANUP] Removing intermediate files..."

rm -f \
    "${BASE}.snp.vcf"              "${BASE}.snp.vcf.idx" \
    "${BASE}.indel.vcf"            "${BASE}.indel.vcf.idx" \
    "${BASE}.snp.3x.vcf"           "${BASE}.snp.3x.vcf.idx" \
    "${BASE}.snp.3x.nocall.vcf"    "${BASE}.snp.3x.nocall.vcf.idx" \
    "${BASE}.snp.3x.filtered.vcf"  "${BASE}.snp.3x.filtered.vcf.idx" \
    "${BASE}.indel.3x.vcf"         "${BASE}.indel.3x.vcf.idx" \
    "${BASE}.indel.3x.nocall.vcf"  "${BASE}.indel.3x.nocall.vcf.idx" \
    "${BASE}.indel.3x.filtered.vcf"  "${BASE}.indel.3x.filtered.vcf.idx" \
    "${BASE}.snp.10x.vcf"           "${BASE}.snp.10x.vcf.idx" \
    "${BASE}.snp.10x.nocall.vcf"    "${BASE}.snp.10x.nocall.vcf.idx" \
    "${BASE}.snp.10x.filtered.vcf"  "${BASE}.snp.10x.filtered.vcf.idx" \
    "${BASE}.indel.10x.vcf"         "${BASE}.indel.10x.vcf.idx" \
    "${BASE}.indel.10x.nocall.vcf"  "${BASE}.indel.10x.nocall.vcf.idx" \
    "${BASE}.indel.10x.filtered.vcf"  "${BASE}.indel.10x.filtered.vcf.idx"

echo "[CLEANUP] Done."

# -----------------------------------------------------------------------------
# SECTION 9: Summary
# -----------------------------------------------------------------------------
echo ""
echo "================================================================"
echo "[DONE] All results in: $OUTDIR"
echo "================================================================"
echo ""
echo "Per-variant CSVs:"
ls -lh "$OUTDIR"/*.csv 2>/dev/null | grep -v count_summary \
    || echo "  (none found)"
echo ""
echo "Count summary: $OUTDIR/ancestor_count_summary.csv"
