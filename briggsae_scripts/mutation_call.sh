#!/usr/bin/env bash

#SBATCH --job-name=mutation_call
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=30gb
#SBATCH --time=24:00:00
#SBATCH --output=mutation_call_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou

set -euo pipefail

# =============================================================================
# Short-read mutation calling + ANC-only site detection
#
# Input  : ./filtered/*.af_filtered.vcf
#          (het sites and >10 missing already removed by AF filter)
#
# Steps per file:
#   1. Call mutations: sites where exactly 1 sample is homozygous alt (1/1)
#   2. Identify ANC-only sites: the single variant sample contains "ANC"
#      → collected across all 8 files into one combined CSV
#   3. Save clean mutation VCF: ANC-only sites removed
#
# Outputs:
#   ./mutations/
#       B1.bp.GVCFs.snp.3x.mutation.vcf.gz         ← all single-variant sites
#       B1.bp.GVCFs.snp.3x.mutation.clean.vcf.gz   ← ANC-only sites removed
#       ... (same for all 8 files)
#   ./mutations/ANC_only_sites.csv  ← combined list from all 8 files
# =============================================================================

module purge
module load bcftools

INPUT_DIR="./filtered"
OUT_DIR="./mutations"
ANC_CSV="${OUT_DIR}/ANC_only_sites.csv"

mkdir -p "$OUT_DIR"

# Write CSV header
echo "CHROM,POS,REF,ALT,ANC_sample,coverage" > "$ANC_CSV"

# =============================================================================
# Process each af_filtered VCF
# =============================================================================

for vcf in "${INPUT_DIR}"/*.af_filtered.vcf; do

    base=$(basename "$vcf")
    echo ""
    echo "============================================================"
    echo "[INFO] Processing: ${base}"
    echo "============================================================"

    # -------------------------------------------------------------------------
    # Extract coverage label (3x or 10x) from filename
    # e.g. B1.bp.GVCFs.snp.3x.filtered.biallelic.af_filtered.vcf -> 3x
    # -------------------------------------------------------------------------
    if [[ "$base" == *".5x."* ]]; then
        cov="5x"
    elif [[ "$base" == *".7x."* ]]; then
        cov="7x"
    else
        echo "[WARN] Could not determine coverage from filename: $base — skipping"
        continue
    fi

    # Build output file prefix:
    # B1.bp.GVCFs.snp.3x.filtered.biallelic.af_filtered.vcf
    #   -> B1.bp.GVCFs.snp.3x
    prefix=$(echo "$base" | sed 's/\.filtered\.biallelic\.af_filtered\.vcf//')
    out_all="${OUT_DIR}/${prefix}.mutation.vcf.gz"
    out_clean="${OUT_DIR}/${prefix}.mutation.clean.vcf.gz"

    # -------------------------------------------------------------------------
    # Step 1: Compress + index the input VCF for bcftools
    # -------------------------------------------------------------------------
    gz="${INPUT_DIR}/${base}.gz"
    if [[ ! -f "$gz" ]]; then
        echo "[INFO] Compressing: ${base}"
        bcftools view -Oz -o "$gz" "$vcf"
    fi
    if [[ ! -f "${gz}.csi" ]]; then
        bcftools index -f "$gz"
    fi

    # -------------------------------------------------------------------------
    # Step 2: Call mutations
    #   - Exactly 1 sample is homozygous alt (1/1)
    #   - Hets and >10 missing already handled by AF filter upstream
    # -------------------------------------------------------------------------
    echo "[INFO] Calling mutations (COUNT(GT=alt)=1)..."
    bcftools view -i 'COUNT(GT="alt")=1' "$gz" -Oz -o "$out_all"
    bcftools index -f "$out_all"

    n_all=$(bcftools view -H "$out_all" | wc -l)
    echo "[INFO] Mutation sites (all): ${n_all}"

    # -------------------------------------------------------------------------
    # Step 3: Get sample names from VCF
    # -------------------------------------------------------------------------
    mapfile -t all_samples < <(bcftools query -l "$gz")

    # -------------------------------------------------------------------------
    # Step 4: Identify ANC-only sites
    #   A site is ANC-only if the single variant sample has "ANC" in its name.
    #   We check this by querying per-sample genotypes and finding which sample
    #   is 1/1 or 1|1 at each site.
    # -------------------------------------------------------------------------
    echo "[INFO] Detecting ANC-only sites..."

    # Build sample list as awk array for name lookup
    # Output: CHROM POS REF ALT [SAMPLE=GT ...]
    # Then awk finds the single 1/1 carrier and checks for ANC in name
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' \
        "$out_all" \
    | awk -v OFS="," -v cov="$cov" '
        function is_homalt(gt) {
            return (gt == "1/1" || gt == "1|1")
        }
        {
            chrom=$1; pos=$2; ref=$3; alt=$4
            carrier=""
            for (i=5; i<=NF; i++) {
                split($i, a, "=")
                samp=a[1]; gt=a[2]
                if (is_homalt(gt)) {
                    carrier=samp
                    break
                }
            }
            # Only flag if carrier name contains "ANC"
            if (carrier != "" && index(carrier, "ANC") > 0) {
                print chrom, pos, ref, alt, carrier, cov
            }
        }
    ' >> "$ANC_CSV"

    # Count how many ANC-only sites were found for this file
    # (count lines added since last iteration — approximate via grep)
    n_anc=$(awk -v OFS="," -v cov="$cov" -v base="$base" '
        NR>1 && $6==cov { count++ }
        END { print count+0 }
    ' "$ANC_CSV")
    echo "[INFO] ANC-only sites detected: ${n_anc}"

    # -------------------------------------------------------------------------
    # Step 5: Build exclusion list of ANC-only positions for this file
    #   then subtract them from the mutation VCF to produce the clean file
    # -------------------------------------------------------------------------

    # Extract ANC-only CHROM+POS for this file from the CSV
    # We re-run the query rather than parsing the CSV to avoid complexity
    anc_sites_tmp="${OUT_DIR}/tmp_anc_sites_${prefix}.txt"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' \
        "$out_all" \
    | awk '
        function is_homalt(gt) { return (gt == "1/1" || gt == "1|1") }
        {
            for (i=5; i<=NF; i++) {
                split($i, a, "=")
                if (is_homalt(a[2]) && index(a[1], "ANC") > 0) {
                    print $1"\t"$2
                    break
                }
            }
        }
    ' > "$anc_sites_tmp"

    n_excl=$(wc -l < "$anc_sites_tmp")
    echo "[INFO] Sites to exclude: ${n_excl}"

    if [[ "$n_excl" -gt 0 ]]; then
        # Use bcftools view -T ^file to exclude those positions
        bcftools view -T "^${anc_sites_tmp}" "$out_all" -Oz -o "$out_clean"
    else
        # No ANC-only sites — clean file is identical to mutation file
        cp "$out_all" "$out_clean"
    fi

    bcftools index -f "$out_clean"
    rm -f "$anc_sites_tmp"

    n_clean=$(bcftools view -H "$out_clean" | wc -l)
    echo "[INFO] Mutation sites after ANC removal: ${n_clean}"
    echo "[INFO] Written: $(basename $out_all)"
    echo "[INFO] Written: $(basename $out_clean)"

done

# =============================================================================
# Final summary
# =============================================================================

echo ""
echo "============================================================"
echo "[DONE] All files processed."
echo ""
echo "ANC-only sites combined: ${ANC_CSV}"
echo "Total ANC-only entries  : $(tail -n +2 "$ANC_CSV" | wc -l)"
echo ""
echo "Output VCFs:"
ls -lh "${OUT_DIR}"/*.vcf.gz 2>/dev/null || echo "  (none found)"
echo "============================================================"
