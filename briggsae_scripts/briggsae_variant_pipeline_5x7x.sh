#!/bin/bash
#SBATCH --job-name=briggsae_s11_s12_5x7x
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=15gb
#SBATCH --time=480:00:00
#SBATCH --output=briggsae_s11_s12_5x7x_%j.out
#SBATCH --account=baer --qos=baer

set -euo pipefail

# =============================================================================
# SETTINGS — edit here only
# =============================================================================
WORKDIR="/orange/baer/briggsae"
REF="20250626_c_briggsae_Feb2020.genome.fa"
BASES=("B1.bp.GVCFs" "B2.bp.GVCFs")   # loop over both datasets

cd "$WORKDIR"
module load gatk/
module load bcftools/

JAVA_MEM="-Xmx8g"

# =============================================================================
# HELPER: run_dp_filter
# NOTE: writes output path to a temp file to avoid capturing GATK stdout
# =============================================================================
run_dp_filter() {
    local INPUT="$1"
    local DP="$2"
    local TAG="$3"
    local RETFILE="$4"
    local BASE="$5"

    local FILTERED="${BASE}.${TAG}.vcf"
    local NOCALL="${BASE}.${TAG}.nocall.vcf"
    local CLEAN="${BASE}.${TAG}.filtered.vcf"

    echo "[STEP 11] DP${DP} filter on $INPUT -> $CLEAN"

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

# =============================================================================
# MAIN LOOP: process B1 and B2
# Steps 9 and 10 are skipped — snp.vcf and indel.vcf already exist
# Only 5x and 7x thresholds are generated here (3x and 10x already done)
# =============================================================================
for BASE in "${BASES[@]}"; do
    echo ""
    echo "============================================================"
    echo " Processing: $BASE"
    echo "============================================================"

    SNP="${BASE}.snp.vcf"
    INDEL="${BASE}.indel.vcf"

    [[ -f "$SNP"   ]] || { echo "ERROR: SNP VCF not found: $SNP"     >&2; exit 1; }
    [[ -f "$INDEL" ]] || { echo "ERROR: Indel VCF not found: $INDEL" >&2; exit 1; }
    echo "[INFO] Using existing: $SNP, $INDEL — skipping Steps 9 & 10."

    # -------------------------------------------------------------------------
    # STEP 11: DP filtering — 5x and 7x only
    # -------------------------------------------------------------------------
    RETFILE_SNP5=".${BASE}_snp5x_out"
    RETFILE_IND5=".${BASE}_ind5x_out"
    RETFILE_SNP7=".${BASE}_snp7x_out"
    RETFILE_IND7=".${BASE}_ind7x_out"

    run_dp_filter "$SNP"   5 "snp.5x"   "$RETFILE_SNP5" "$BASE"
    run_dp_filter "$INDEL" 5 "indel.5x" "$RETFILE_IND5" "$BASE"
    run_dp_filter "$SNP"   7 "snp.7x"   "$RETFILE_SNP7" "$BASE"
    run_dp_filter "$INDEL" 7 "indel.7x" "$RETFILE_IND7" "$BASE"

    SNP_5X_FILTERED=$(cat "$RETFILE_SNP5")
    INDEL_5X_FILTERED=$(cat "$RETFILE_IND5")
    SNP_7X_FILTERED=$(cat "$RETFILE_SNP7")
    INDEL_7X_FILTERED=$(cat "$RETFILE_IND7")
    rm -f "$RETFILE_SNP5" "$RETFILE_IND5" "$RETFILE_SNP7" "$RETFILE_IND7"

    echo "[STEP 11] Done."

    # -------------------------------------------------------------------------
    # STEP 12: Biallelic filtering
    # -------------------------------------------------------------------------
    echo "[STEP 12] Applying biallelic filter..."

    SNP_5X_FINAL="${BASE}.snp.5x.filtered.biallelic.vcf"
    INDEL_5X_FINAL="${BASE}.indel.5x.filtered.biallelic.vcf"
    SNP_7X_FINAL="${BASE}.snp.7x.filtered.biallelic.vcf"
    INDEL_7X_FINAL="${BASE}.indel.7x.filtered.biallelic.vcf"

    bcftools view --max-alleles 2 "$SNP_5X_FILTERED"   -o "$SNP_5X_FINAL"
    bcftools view --max-alleles 2 "$INDEL_5X_FILTERED" -o "$INDEL_5X_FINAL"
    bcftools view --max-alleles 2 "$SNP_7X_FILTERED"   -o "$SNP_7X_FINAL"
    bcftools view --max-alleles 2 "$INDEL_7X_FILTERED" -o "$INDEL_7X_FINAL"

    echo "[STEP 12] Done."

    # -------------------------------------------------------------------------
    # CLEANUP: intermediates for this BASE only
    # -------------------------------------------------------------------------
    echo "[CLEANUP] Removing intermediates for $BASE..."
    rm -f \
        "${BASE}.snp.5x.vcf"            "${BASE}.snp.5x.vcf.idx" \
        "${BASE}.snp.5x.nocall.vcf"     "${BASE}.snp.5x.nocall.vcf.idx" \
        "${BASE}.snp.5x.filtered.vcf"   "${BASE}.snp.5x.filtered.vcf.idx" \
        "${BASE}.indel.5x.vcf"          "${BASE}.indel.5x.vcf.idx" \
        "${BASE}.indel.5x.nocall.vcf"   "${BASE}.indel.5x.nocall.vcf.idx" \
        "${BASE}.indel.5x.filtered.vcf" "${BASE}.indel.5x.filtered.vcf.idx" \
        "${BASE}.snp.7x.vcf"            "${BASE}.snp.7x.vcf.idx" \
        "${BASE}.snp.7x.nocall.vcf"     "${BASE}.snp.7x.nocall.vcf.idx" \
        "${BASE}.snp.7x.filtered.vcf"   "${BASE}.snp.7x.filtered.vcf.idx" \
        "${BASE}.indel.7x.vcf"          "${BASE}.indel.7x.vcf.idx" \
        "${BASE}.indel.7x.nocall.vcf"   "${BASE}.indel.7x.nocall.vcf.idx" \
        "${BASE}.indel.7x.filtered.vcf" "${BASE}.indel.7x.filtered.vcf.idx"
    echo "[CLEANUP] Done."

    # -------------------------------------------------------------------------
    # SUMMARY for this BASE
    # -------------------------------------------------------------------------
    echo ""
    echo " Files produced for $BASE:"
    echo "  $SNP_5X_FINAL"
    echo "  $INDEL_5X_FINAL"
    echo "  $SNP_7X_FINAL"
    echo "  $INDEL_7X_FINAL"

done

echo ""
echo "============================================================"
echo " All done."
echo "============================================================"
