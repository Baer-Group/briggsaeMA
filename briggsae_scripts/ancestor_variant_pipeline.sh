#!/usr/bin/env bash
#SBATCH --job-name=anc_variants
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --time=96:00:00
#SBATCH --output=anc_variants_%j.out
#SBATCH --account=baer --qos=baer

set -euo pipefail

# =============================================================================
# SETTINGS — edit here only
# =============================================================================
WORKDIR="/blue/baer/m.rifat/briggsae_ancestors"
REF="20250626_c_briggsae_Feb2020.genome.fa"
BASE="briggsae_ancestor.bp.GVCFs"   # base name for all output files

cd "$WORKDIR"
module load gatk/
module load bcftools/
module load python/3.10

JAVA_MEM="-Xmx39g"

# This script assumes briggsae_ancestor.bp.GVCFs.vcf already exists.
JOINT="${BASE}.vcf"
[[ -f "$JOINT" ]] || { echo "ERROR: Joint VCF not found: $JOINT" >&2; exit 1; }
echo "[INFO] Joint VCF found: $JOINT"

# =============================================================================
# STEP 1: Split SNPs and indels
# =============================================================================
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

# =============================================================================
# STEP 2: DP filtering (3x and 10x) for each variant type
# NOTE: run_dp_filter writes its output path to a temp file rather than
# using $() command substitution to avoid capturing GATK verbose logs.
# =============================================================================

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

# =============================================================================
# STEP 3: Biallelic filter
# =============================================================================
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

# =============================================================================
# STEP 4: Per-file Python filtering + CSV export
#
#   For each biallelic VCF:
#     - For each 1/1 sample: check minority AF (REF_depth/total_depth)
#       If minority AF > 10% → mask that sample to ./.
#     - Heterozygotes (0/1) are kept as-is, no AF check
#     - Drop site if no 1/1 remains after masking
#     - Drop site if more than 1 sample is missing after masking
#     - Output: CSV with CHROM, POS, REF, ALT, [sample_name cols]
# =============================================================================
echo "[STEP 4] Applying AF + missing filters and exporting CSVs..."

process_vcf() {
    local INPUT_VCF="$1"
    local OUTPUT_CSV="$2"

python3 << PYEOF
import sys

input_vcf  = "${INPUT_VCF}"
output_csv = "${OUTPUT_CSV}"

MAX_MISSING  = 1      # at most 1 missing genotype allowed across 6 ancestors
AF_THRESHOLD = 0.10   # minority AF threshold for 1/1 samples

sites_total    = 0
sites_no_hom   = 0
sites_missing  = 0
sites_kept     = 0

samples = []
header_written = False

with open(input_vcf, 'r') as fin, open(output_csv, 'w') as fout:
    for line in fin:
        # ------------------------------------------------------------------
        # Parse header to get sample names
        # ------------------------------------------------------------------
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            fields = line.rstrip('\n').split('\t')
            samples = fields[9:]
            # Write CSV header
            fout.write('CHROM,POS,REF,ALT,' + ','.join(samples) + '\n')
            header_written = True
            continue

        if not header_written:
            continue

        sites_total += 1
        fields = line.rstrip('\n').split('\t')
        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]

        format_fields = fields[8].split(':')
        try:
            gt_idx = format_fields.index('GT')
        except ValueError:
            sites_total -= 1
            continue
        try:
            ad_idx = format_fields.index('AD')
        except ValueError:
            ad_idx = None

        sample_data = fields[9:]

        # ------------------------------------------------------------------
        # Per-sample processing:
        #   - Only 1/1 samples get AF check
        #   - 0/1 kept as-is
        #   - 0/0 kept as-is
        #   - ./. kept as-is
        # ------------------------------------------------------------------
        new_samples = []
        for s in sample_data:
            sf = s.split(':')
            if gt_idx >= len(sf):
                new_samples.append(s)
                continue

            gt = sf[gt_idx].replace('|', '/')
            alleles = gt.split('/')

            # Only apply AF masking to homozygous alt (1/1)
            if set(alleles) == {'1'} and ad_idx is not None:
                if ad_idx < len(sf) and sf[ad_idx] not in ('.', ''):
                    try:
                        ad_vals = [int(x) for x in sf[ad_idx].split(',')]
                        total_depth = sum(ad_vals)
                        if total_depth > 0:
                            ref_depth = ad_vals[0]
                            minority_af = ref_depth / total_depth
                            if minority_af > AF_THRESHOLD:
                                sf[gt_idx] = './.'
                                new_samples.append(':'.join(sf))
                                continue
                    except ValueError:
                        pass

            new_samples.append(s)

        # ------------------------------------------------------------------
        # Filter 1: Must have at least one 1/1 remaining after AF masking
        # ------------------------------------------------------------------
        has_hom_alt = False
        for s in new_samples:
            sf = s.split(':')
            if gt_idx >= len(sf):
                continue
            gt = sf[gt_idx].replace('|', '/')
            if set(gt.split('/')) == {'1'}:
                has_hom_alt = True
                break

        if not has_hom_alt:
            sites_no_hom += 1
            continue

        # ------------------------------------------------------------------
        # Filter 2: At most 1 missing genotype allowed
        # ------------------------------------------------------------------
        missing_count = sum(
            1 for s in new_samples
            if '.' in s.split(':')[gt_idx].replace('|', '/').split('/')
            if gt_idx < len(s.split(':'))
        )

        if missing_count > MAX_MISSING:
            sites_missing += 1
            continue

        # ------------------------------------------------------------------
        # Site passes — write to CSV
        # Extract GT only for the CSV columns
        # ------------------------------------------------------------------
        gt_vals = []
        for s in new_samples:
            sf = s.split(':')
            gt_vals.append(sf[gt_idx] if gt_idx < len(sf) else './.')

        fout.write(','.join([chrom, pos, ref, alt] + gt_vals) + '\n')
        sites_kept += 1

print(f"  Input VCF      : {input_vcf}")
print(f"  Total sites    : {sites_total}")
print(f"  No hom-alt     : {sites_no_hom}")
print(f"  >1 missing     : {sites_missing}")
print(f"  Sites kept     : {sites_kept}")
print(f"  Output CSV     : {output_csv}")
PYEOF
}

process_vcf "$SNP_3X_BI"    "${BASE}.snp.3x.csv"
process_vcf "$INDEL_3X_BI"  "${BASE}.indel.3x.csv"
process_vcf "$SNP_10X_BI"   "${BASE}.snp.10x.csv"
process_vcf "$INDEL_10X_BI" "${BASE}.indel.10x.csv"

echo "[STEP 4] Done."

# =============================================================================
# CLEANUP: Remove intermediates, keep joint VCF, biallelic VCFs, and CSVs
# =============================================================================
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

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo " Pipeline complete. Files retained:"
echo "  [joint]     $JOINT"
echo "  [biallelic] $SNP_3X_BI"
echo "  [biallelic] $INDEL_3X_BI"
echo "  [biallelic] $SNP_10X_BI"
echo "  [biallelic] $INDEL_10X_BI"
echo "  [csv]       ${BASE}.snp.3x.csv"
echo "  [csv]       ${BASE}.indel.3x.csv"
echo "  [csv]       ${BASE}.snp.10x.csv"
echo "  [csv]       ${BASE}.indel.10x.csv"
echo "============================================================"
