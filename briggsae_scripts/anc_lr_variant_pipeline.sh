#!/usr/bin/env bash
#SBATCH --job-name=anc_lr_variants
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --time=120:00:00
#SBATCH --output=anc_lr_variants_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou

# =============================================================================
# C. briggsae Long-Read Ancestor Variant Listing Pipeline
# =============================================================================
# OVERVIEW:
#   1. Index per-sample gVCF files (naming: *.dv.g.renamed.vcf.gz)
#   2. Joint-call using GLnexus
#   3. Convert BCF → VCF
#   4. For each variant type (SNP / INDEL) and coverage threshold (3x / 10x):
#       a. Extract variant type
#       b. Mark low-DP genotypes as missing
#       c. Keep only biallelic sites
#       d. AF filter: for 1/1 samples only, mask to ./. if minority AF > 10%
#          Heterozygotes are kept as-is. No het-site removal.
#       e. Keep sites with at least one 1/1 remaining after AF masking
#          (with 2 samples, missing filter is redundant — if one sample is
#           missing and the other is not 1/1, the site is dropped anyway)
#       f. Export to CSV: CHROM, POS, REF, ALT, [sample_name columns]
#
# NOTE ON MODULE LOADING:
#   glnexus and python/3.10 conflict on HiPerGator — loaded/unloaded at point
#   of use to avoid module replacement.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# SECTION 1: User settings
# -----------------------------------------------------------------------------

WORKDIR="/blue/baer/m.rifat/briggsae_ancestors"
REF="20250626_c_briggsae_Feb2020.genome.fa"
OUTDIR="anc_lr_out"
THREADS=8

DP_LEVELS=(3 10)

AF_THRESHOLD=0.10

# -----------------------------------------------------------------------------
# SECTION 2: Initial setup
# -----------------------------------------------------------------------------

module purge || true
module load bcftools
module load samtools

cd "$WORKDIR"
mkdir -p "$OUTDIR"

[[ -f "$REF" ]] || { echo "ERROR: Reference not found: $REF" >&2; exit 1; }
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

# -----------------------------------------------------------------------------
# SECTION 3: Build gVCF list and index
# -----------------------------------------------------------------------------

echo "[INFO] Building gVCF list..."

GVCF_LIST="gvcfs.list"
ls *.dv.g.renamed.vcf.gz > "$GVCF_LIST"

while read -r gvcf; do
    [[ -f "${gvcf}.tbi" ]] || bcftools index -t "$gvcf"
done < "$GVCF_LIST"

echo "[INFO] Found $(wc -l < "$GVCF_LIST") gVCF files."

# -----------------------------------------------------------------------------
# SECTION 4: Joint variant calling with GLnexus
# -----------------------------------------------------------------------------

module load glnexus

BCF="cohort.anc.bcf"

if [[ -f "$BCF" && -f "${BCF}.csi" ]]; then
    echo "[INFO] BCF already exists, skipping GLnexus: $BCF"
else
    echo "[INFO] Running GLnexus joint calling..."
    rm -rf GLnexus.DB
    glnexus_cli --config DeepVariantWGS \
                --list "$GVCF_LIST" \
        | bcftools view -Ob -o "$BCF"
    bcftools index "$BCF"
    echo "[INFO] Joint calling complete: $BCF"
fi

module unload glnexus

# -----------------------------------------------------------------------------
# SECTION 5: Convert BCF → compressed VCF
# -----------------------------------------------------------------------------

JOINT_VCF="cohort.anc.vcf.gz"

bcftools view "$BCF" -Oz -o "$JOINT_VCF"
bcftools index -t "$JOINT_VCF"
echo "[INFO] Joint VCF ready: $JOINT_VCF"

# -----------------------------------------------------------------------------
# SECTION 6: AF filter function
# -----------------------------------------------------------------------------
# Applies minority AF masking to 1/1 samples only.
# Heterozygotes are passed through unchanged.
# 0/0 samples are passed through unchanged.

run_af_filter() {
    local in_gz="$1"
    local out_gz="$2"
    local af_thresh="$3"
    local tag="$4"

    module load python/3.10

python3 << PYEOF
import gzip

in_path   = "${in_gz}"
out_path  = "${out_gz}"
af_thresh = float("${af_thresh}")
label     = "${tag}"

sites_total    = 0
gt_set_missing = 0
sites_written  = 0

with gzip.open(in_path, 'rt') as fin, gzip.open(out_path, 'wt') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue

        sites_total += 1
        fields = line.rstrip('\n').split('\t')
        fmt    = fields[8].split(':')

        try:
            gt_idx = fmt.index('GT')
            ad_idx = fmt.index('AD')
        except ValueError:
            fout.write(line)
            sites_written += 1
            continue

        samples = fields[9:]
        new_samples = []

        for s in samples:
            sf  = s.split(':')
            gt  = sf[gt_idx].replace('|', '/') if gt_idx < len(sf) else './.'
            als = gt.split('/')

            # Only apply AF masking to homozygous alt (1/1)
            # Hets (0/1) and hom-ref (0/0) pass through unchanged
            if set(als) == {'1'}:
                if ad_idx < len(sf) and sf[ad_idx] not in ('.', ''):
                    try:
                        ad_vals = [int(x) for x in sf[ad_idx].split(',')]
                        total   = sum(ad_vals)
                        if total > 0:
                            minority_af = ad_vals[0] / total  # REF depth / total
                            if minority_af > af_thresh:
                                sf[gt_idx] = './.'
                                new_samples.append(':'.join(sf))
                                gt_set_missing += 1
                                continue
                    except ValueError:
                        pass

            new_samples.append(s)

        fields[9:] = new_samples
        fout.write('\t'.join(fields) + '\n')
        sites_written += 1

print(f"  [{label}] AF filter summary:")
print(f"    Sites evaluated          : {sites_total}")
print(f"    Genotypes set to missing : {gt_set_missing}")
print(f"    Sites written            : {sites_written}", flush=True)
PYEOF

    module unload python/3.10
}

# -----------------------------------------------------------------------------
# SECTION 7: CSV export function
# -----------------------------------------------------------------------------
# Reads a VCF (standard gzip), keeps sites with at least one 1/1 genotype,
# outputs CHROM, POS, REF, ALT, [sample GT columns].

export_csv() {
    local in_gz="$1"
    local out_csv="$2"
    local tag="$3"

    module load python/3.10

python3 << PYEOF
import gzip

in_path  = "${in_gz}"
out_path = "${out_csv}"
label    = "${tag}"

sites_total   = 0
sites_no_hom  = 0
sites_kept    = 0

with gzip.open(in_path, 'rt') as fin, open(out_path, 'w') as fout:
    samples = []
    for line in fin:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            samples = line.rstrip('\n').split('\t')[9:]
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

        sample_data = fields[9:]

        # Keep site only if at least one sample has 1/1
        has_hom_alt = any(
            set(s.split(':')[gt_idx].replace('|', '/').split('/')) == {'1'}
            for s in sample_data
            if gt_idx < len(s.split(':'))
        )

        if not has_hom_alt:
            sites_no_hom += 1
            continue

        # Extract GT only for CSV columns
        gt_vals = [
            s.split(':')[gt_idx] if gt_idx < len(s.split(':')) else './.'
            for s in sample_data
        ]

        fout.write(','.join([chrom, pos, ref, alt] + gt_vals) + '\n')
        sites_kept += 1

print(f"  [{label}] CSV export summary:")
print(f"    Sites evaluated : {sites_total}")
print(f"    No hom-alt      : {sites_no_hom}")
print(f"    Sites kept      : {sites_kept}")
print(f"    Output CSV      : {out_path}", flush=True)
PYEOF

    module unload python/3.10
}

# -----------------------------------------------------------------------------
# SECTION 8: Per-type, per-depth processing
# -----------------------------------------------------------------------------

process_variants() {
    local vartype="$1"
    local dp="$2"
    local tag="${vartype}.DP${dp}"

    echo ""
    echo "============================================================"
    echo "[INFO] Processing: ${vartype} | DP >= ${dp}"
    echo "============================================================"

    # Step a-c: extract type → mask low-DP → biallelic
    local tmp_biallelic="$OUTDIR/${tag}.biallelic.vcf.gz"

    bcftools view -v "$vartype" "$JOINT_VCF" \
    | bcftools +setGT -- -t q -n . -i "FMT/DP<${dp}" \
    | bcftools view -m2 -M2 -Oz -o "$tmp_biallelic"

    bcftools index -t "$tmp_biallelic"

    # Step d: AF filter (1/1 only, hets unchanged)
    local tmp_af="$OUTDIR/${tag}.af_filtered.vcf.gz"
    echo "[INFO] Applying AF filter..."
    run_af_filter "$tmp_biallelic" "$tmp_af" "$AF_THRESHOLD" "$tag"

    # Step e+f: keep sites with ≥1 hom-alt, export CSV
    local out_csv="$OUTDIR/${tag}.ancestor_variants.csv"
    echo "[INFO] Exporting CSV..."
    export_csv "$tmp_af" "$out_csv" "$tag"

    echo "[INFO] Done: $out_csv"
}

# -----------------------------------------------------------------------------
# SECTION 9: Run all type × depth combinations
# -----------------------------------------------------------------------------

for dp in "${DP_LEVELS[@]}"; do
    process_variants "snps"   "$dp"
    process_variants "indels" "$dp"
done

# -----------------------------------------------------------------------------
# SECTION 10: Per-sample hom-alt counts across all CSVs
# -----------------------------------------------------------------------------

echo ""
echo "[INFO] Counting hom-alt (1/1) sites per sample..."

module load python/3.10

python3 << PYEOF
import os, csv
from collections import defaultdict

outdir  = "${OUTDIR}"
out_csv = os.path.join(outdir, "per_sample_counts.csv")

files = sorted(
    f for f in os.listdir(outdir)
    if f.endswith('.ancestor_variants.csv')
)

if not files:
    print("  No ancestor_variants.csv files found.")
else:
    all_tags    = []
    all_samples = []
    counts      = {}   # tag -> {sample: int}

    shared = {}   # tag -> count of sites where ALL samples are 1/1

    for fname in files:
        tag = fname.replace('.ancestor_variants.csv', '')
        all_tags.append(tag)
        counts[tag]  = defaultdict(int)
        shared[tag]  = 0

        fpath = os.path.join(outdir, fname)
        with open(fpath, 'r') as f:
            reader = csv.DictReader(f)
            samples = [c for c in reader.fieldnames if c not in ('CHROM','POS','REF','ALT')]
            if not all_samples:
                all_samples = samples
            for row in reader:
                hom_flags = []
                for samp in samples:
                    gt = row[samp].replace('|', '/')
                    is_hom = set(gt.split('/')) == {'1'}
                    if is_hom:
                        counts[tag][samp] += 1
                    hom_flags.append(is_hom)
                # Shared: all samples are 1/1 at this site
                if all(hom_flags):
                    shared[tag] += 1



    print(f"\n  Written: {out_csv}", flush=True)
PYEOF

module unload python/3.10

# -----------------------------------------------------------------------------
# SECTION 11: Summary
# -----------------------------------------------------------------------------

echo ""
echo "================================================================"
echo "[DONE] All results in: $OUTDIR"
echo "================================================================"
echo ""
echo "Output CSVs:"
ls -lh "$OUTDIR"/*.ancestor_variants.csv 2>/dev/null \
    || echo "  (none found — check log for errors)"
echo ""
echo "Per-sample counts: $OUTDIR/per_sample_counts.csv"
