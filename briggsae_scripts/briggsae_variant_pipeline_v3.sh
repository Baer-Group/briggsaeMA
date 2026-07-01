#!/usr/bin/env bash
#SBATCH --job-name=dv_pipeline
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --time=120:00:00
#SBATCH --output=dv_pipeline_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou

# =============================================================================
# C. briggsae DeepVariant Cohort Variant-Calling Pipeline
# =============================================================================
# OVERVIEW OF STEPS:
#   1. Index all per-sample gVCF files
#   2. Joint-call variants across all samples using GLnexus
#   3. Convert the joint-call BCF output to VCF
#   4. For each variant type (SNP / INDEL) and each coverage threshold (3x / 10x):
#       a. Extract variant type
#       b. Mark low-coverage genotypes as missing (instead of dropping the site)
#       c. Keep only strictly biallelic sites (1 REF + 1 ALT)
#          Note: normalization skipped for consistency with short-read pipeline
#       e. AF filter: discard het sites; set minority AF > 10% genotypes to missing
#       f. Drop sites with 3 or more missing genotypes (allow up to 2)
#       g. Keep mutation sites: exactly 1 homozygous-alt sample, 0 heterozygotes
#       h. Generate bcftools stats report
#
# NOTE ON MODULE LOADING:
#   glnexus and python/3.10 conflict on HiPerGator (Lmod replaces one with the
#   other if loaded together). To avoid this, modules are loaded and unloaded
#   at the point they are needed rather than all upfront:
#     Step 2  : glnexus only
#     Step 4e : python/3.10 only  (AF filter)
#     All else: bcftools + samtools only
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# SECTION 1: User settings
# -----------------------------------------------------------------------------

WORKDIR="/orange/juannanzhou/Rifat_CB_raw/briggsae/Long_read_Fresh"
REF="20250626_c_briggsae_Feb2020.genome.fa"
OUTDIR="dv_pipeline_out_v2"
THREADS=8

#DP_LEVELS=(3 10)
DP_LEVELS=(5 7)

MUTATION_EXPR='COUNT(GT="alt")=1 && COUNT(GT="het")=0'

AF_THRESHOLD=0.10
MAX_MISSING=3

# -----------------------------------------------------------------------------
# SECTION 2: Initial setup — load only what is needed for setup steps
# -----------------------------------------------------------------------------

module purge || true
module load bcftools
module load samtools

cd "$WORKDIR"
mkdir -p "$OUTDIR"

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found: $REF" >&2
    exit 1
fi

if [[ ! -f "${REF}.fai" ]]; then
    echo "[INFO] Indexing reference FASTA..."
    samtools faidx "$REF"
fi

# -----------------------------------------------------------------------------
# SECTION 3: Build the gVCF list and index each gVCF
# -----------------------------------------------------------------------------

echo "[INFO] Building gVCF list and checking indexes..."

GVCF_LIST="gvcfs.list"
ls *.dv.g.renamed.vcf.gz > "$GVCF_LIST"

while read -r gvcf; do
    if [[ ! -f "${gvcf}.tbi" ]]; then
        echo "[INFO] Indexing ${gvcf}..."
        bcftools index -t "$gvcf"
    fi
done < "$GVCF_LIST"

echo "[INFO] Found $(wc -l < "$GVCF_LIST") gVCF files."

# -----------------------------------------------------------------------------
# SECTION 4: Joint variant calling with GLnexus
# -----------------------------------------------------------------------------
# glnexus conflicts with python/3.10 — load it here, unload before Python step

echo "[INFO] Loading glnexus for joint calling..."
module load glnexus

echo "[INFO] Running GLnexus joint calling..."

BCF="cohort.dvglnexus.bcf"

glnexus_cli --config DeepVariantWGS \
            --list "$GVCF_LIST" \
    | bcftools view -Ob -o "$BCF"

bcftools index "$BCF"

echo "[INFO] Joint calling complete: $BCF"

# Unload glnexus immediately after use to free the module slot
module unload glnexus
echo "[INFO] glnexus unloaded."

# -----------------------------------------------------------------------------
# SECTION 5: Convert BCF to compressed VCF
# -----------------------------------------------------------------------------

echo "[INFO] Converting BCF to VCF..."

JOINT_VCF="cohort.dvglnexus.vcf.gz"

bcftools view "$BCF" -Oz -o "$JOINT_VCF"
bcftools index -t "$JOINT_VCF"

echo "[INFO] Joint VCF ready: $JOINT_VCF"

# -----------------------------------------------------------------------------
# SECTION 6: AF filter function
# -----------------------------------------------------------------------------
# Called once per tag. Loads python/3.10, runs the filter, then unloads it
# so bcftools remains the active module for subsequent steps.
#
# Filters applied:
#   1. Het site filter  : discard entire site if any sample is heterozygous
#   2. Minority AF check: set genotype to ./. if minority AF > AF_THRESHOLD
#        - 0/0 call: minority AF = ALT_depth / total
#        - 1/1 call: minority AF = REF_depth / total

run_af_filter() {
    local in_gz="$1"
    local out_gz="$2"
    local af_thresh="$3"
    local tag="$4"

    echo "[INFO] Loading python/3.10 for AF filter..."
    module load python/3.10

    python3 << PYEOF
import gzip

in_path   = "${in_gz}"
out_path  = "${out_gz}"
af_thresh = float("${af_thresh}")
label     = "${tag}"

sites_total    = 0
sites_het      = 0
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

        # Filter 1: discard site if any sample is heterozygous
        het_found = False
        for s in samples:
            sf  = s.split(':')
            gt  = sf[gt_idx].replace('|', '/') if gt_idx < len(sf) else './.'
            als = gt.split('/')
            if '.' not in als and len(set(als)) > 1:
                het_found = True
                break
        if het_found:
            sites_het += 1
            continue

        # Filter 2: per-sample minority AF check
        new_samples = []
        for s in samples:
            sf = s.split(':')
            gt = sf[gt_idx].replace('|', '/') if gt_idx < len(sf) else './.'
            als = gt.split('/')

            if '.' in als:
                new_samples.append(s)
                continue

            if ad_idx >= len(sf) or sf[ad_idx] in ('.', ''):
                new_samples.append(s)
                continue

            try:
                ad_vals = [int(x) for x in sf[ad_idx].split(',')]
            except ValueError:
                new_samples.append(s)
                continue

            total = sum(ad_vals)
            if total == 0:
                sf[gt_idx] = './.'
                new_samples.append(':'.join(sf))
                gt_set_missing += 1
                continue

            allele_set = set(als)

            if allele_set == {'0'}:
                minority_af = sum(ad_vals[1:]) / total if len(ad_vals) > 1 else 0.0
            elif allele_set == {'1'}:
                minority_af = ad_vals[0] / total
            else:
                new_samples.append(s)
                continue

            if minority_af > af_thresh:
                sf[gt_idx] = './.'
                new_samples.append(':'.join(sf))
                gt_set_missing += 1
            else:
                new_samples.append(s)

        fields[9:] = new_samples
        fout.write('\t'.join(fields) + '\n')
        sites_written += 1

print(f"  [{label}] AF filter summary:")
print(f"    Sites evaluated          : {sites_total}")
print(f"    Sites removed (het)      : {sites_het}")
print(f"    Genotypes set to missing : {gt_set_missing}")
print(f"    Sites written            : {sites_written}", flush=True)
PYEOF

    echo "[INFO] Unloading python/3.10..."
    module unload python/3.10
}

# -----------------------------------------------------------------------------
# SECTION 7: Per-type, per-depth processing function
# -----------------------------------------------------------------------------

process_variants() {

    local vartype="$1"
    local dp="$2"
    local tag="${vartype}.DP${dp}"

    echo ""
    echo "============================================================"
    echo "[INFO] Processing: ${vartype} | min coverage DP >= ${dp}"
    echo "============================================================"

    # Steps 4a-4c: extract type → mask low-DP → biallelic
    # Note: normalization (bcftools norm) is intentionally omitted to keep
    # this pipeline consistent with the short-read pipeline, which also uses
    # biallelic filtering only. Multiallelic sites are dropped by -m2 -M2.
    local tmp_biallelic="$OUTDIR/${tag}.biallelic.vcf.gz"

    bcftools view -v "$vartype" "$JOINT_VCF" \
    | bcftools +setGT -- -t q -n . -i "FMT/DP<${dp}" \
    | bcftools view -m2 -M2 -Oz -o "$tmp_biallelic"

    bcftools index -t "$tmp_biallelic"

    # Step 4e: AF filter (loads/unloads python/3.10 internally)
    local tmp_af="$OUTDIR/${tag}.af_filtered.vcf.gz"

    echo "[INFO] Applying AF filter (threshold=${AF_THRESHOLD})..."
    run_af_filter "$tmp_biallelic" "$tmp_af" "$AF_THRESHOLD" "$tag"
    # Note: no index needed here — downstream steps read sequentially.
    # Python's gzip produces standard gzip (not BGZF), which bcftools
    # can read sequentially but cannot index.

    # Step 4f: Remove sites with >= MAX_MISSING missing genotypes
    local tmp_nomissing="$OUTDIR/${tag}.max${MAX_MISSING}missing.vcf.gz"

    echo "[INFO] Filtering sites with >= ${MAX_MISSING} missing genotypes..."
    bcftools view -i "COUNT(GT=\"mis\")<${MAX_MISSING}" \
        "$tmp_af" -Oz -o "$tmp_nomissing"
    bcftools index -t "$tmp_nomissing"

    local n_before n_after
    n_before=$(bcftools view -H "$tmp_af"        | wc -l)
    n_after=$( bcftools view -H "$tmp_nomissing" | wc -l)
    echo "[INFO] Sites before missing filter : ${n_before}"
    echo "[INFO] Sites after missing filter  : ${n_after}"

    # Step 4g: Retain only mutation sites
    local out="$OUTDIR/${tag}.mutations.vcf.gz"

    bcftools view -i "$MUTATION_EXPR" "$tmp_nomissing" -Oz -o "$out"
    bcftools index -t "$out"

    # Step 4h: Generate per-sample statistics
    local stats="$OUTDIR/${tag}.stats.txt"
    bcftools stats -s - "$out" > "$stats"

    local n_sites
    n_sites=$(bcftools view -H "$out" | wc -l)
    echo "[INFO] Output VCF : $out"
    echo "[INFO] Stats file : $stats"
    echo "[INFO] Mutation sites retained: ${n_sites}"
}

# -----------------------------------------------------------------------------
# SECTION 8: Run for all type × depth combinations
# -----------------------------------------------------------------------------

for dp in "${DP_LEVELS[@]}"; do
    process_variants "snps"   "$dp"
    process_variants "indels" "$dp"
done

# -----------------------------------------------------------------------------
# SECTION 9: Final summary
# -----------------------------------------------------------------------------

echo ""
echo "================================================================"
echo "[DONE] Pipeline complete. All results in: $OUTDIR"
echo "================================================================"
echo ""
echo "Output files:"
ls -lh "$OUTDIR"/*.mutations.vcf.gz 2>/dev/null \
    || echo "  (none found — check log for errors)"

echo ""
echo "Quick per-sample mutation counts (from PSC lines in stats files):"
grep '^PSC' "$OUTDIR"/*.stats.txt \
    | awk '{print $1, $4, "hom-alt:", $14}' \
    | column -t
