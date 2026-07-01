#!/usr/bin/env bash
#SBATCH --job-name=contam_matrix
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=02:00:00
#SBATCH --output=contam_matrix_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou

# =============================================================================
# Contamination site detection and variant-sharing matrix
# C. briggsae long-read DeepVariant cohort
# =============================================================================
#
# PURPOSE:
#   After joint variant calling, some sites will show homozygous-alt (1/1)
#   genotypes in MORE THAN ONE line. These multi-line shared sites are likely:
#     - Contamination events (same variant appearing in multiple samples)
#     - Ancestral polymorphisms not fixed in the ancestor
#     - Systematic sequencing/mapping errors
#   This script identifies those sites and summarises how often each pair of
#   lines shares a variant — the off-diagonal values in the matrix are the
#   key signal for contamination between specific sample pairs.
#
# INPUTS:
#   dv_pipeline_out/snps.DP3.max3missing.vcf.gz     (and DP10 equivalents)
#   dv_pipeline_out/indels.DP3.max3missing.vcf.gz   (and DP10 equivalents)
#   These are the direct outputs of briggsae_variant_pipeline_v3.sh, after
#   the AF filter and missing genotype filter have been applied.
#
#   NOTE: If you later decide to run the GQ20 filter (gq20_filter_all.sh),
#   simply update INDIR to "dv_pipeline_out_GQ20" and FILE_SUFFIX to
#   "max3missing.GQ20.vcf.gz" at the top of this script.
#
# OUTPUTS (all written to a single folder: contam_out/):
#   contam_out/contaminant_sites_snps_DP3.csv
#   contam_out/variant_line_matrix_snps_DP3.csv
#   (same pattern for snps_DP10, indels_DP3, indels_DP10)
#
# ANCESTOR FILTER RULE:
#   A site is only considered if BOTH ancestors are homozygous-ref (0/0 or 0|0)
#   with genotype quality GQ >= 20. Sites where either ancestor carries the
#   ALT allele are dropped — these are likely standing polymorphisms in the
#   ancestral background, not new mutations.
#
# CONTAMINATION RULE:
#   A site is flagged as a contaminant/shared site if >= 2 non-ancestor lines
#   are homozygous-alt (1/1 or 1|1) at that site.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# SECTION 1: User settings
# -----------------------------------------------------------------------------

WORKDIR="/orange/juannanzhou/Rifat_CB_raw/briggsae/Long_read_Fresh"
INDIR="dv_pipeline_out"        # direct output of briggsae_variant_pipeline_v3.sh
FILE_SUFFIX="max3missing.vcf.gz"   # update to "max3missing.GQ20.vcf.gz" if GQ20 filter was run
OUTDIR="contam_out"            # single output folder for all results

ANC1="PB_G0_2"
ANC2="HK_G0_1"

ANC_MIN_GQ=20
SEP=";"

# -----------------------------------------------------------------------------
# SECTION 2: Load modules
# -----------------------------------------------------------------------------

module purge || true
module load bcftools || true

# -----------------------------------------------------------------------------
# SECTION 3: Setup
# -----------------------------------------------------------------------------

cd "$WORKDIR"
mkdir -p "$OUTDIR"

# Sanity check — input files must exist
for tag in snps.DP3 snps.DP10 indels.DP3 indels.DP10; do
    f="${INDIR}/${tag}.${FILE_SUFFIX}"
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Input file not found: $f" >&2
        echo "       Run briggsae_variant_pipeline_v3.sh first." >&2
        exit 1
    fi
done

echo "[INFO] Working directory : $WORKDIR"
echo "[INFO] Input directory   : $INDIR"
echo "[INFO] File suffix       : $FILE_SUFFIX"
echo "[INFO] Output directory  : $OUTDIR"
echo "[INFO] Ancestor 1        : $ANC1"
echo "[INFO] Ancestor 2        : $ANC2"
echo "[INFO] Ancestor min GQ   : $ANC_MIN_GQ"

# -----------------------------------------------------------------------------
# SECTION 4: Core processing function
# -----------------------------------------------------------------------------

process_one() {

    local tag="$1"
    local vcf="${INDIR}/${tag}.${FILE_SUFFIX}"

    # Convert dots to underscores for cleaner filenames: snps.DP3 → snps_DP3
    local filetag="${tag//./_}"
    local contam_csv="${OUTDIR}/contaminant_sites_${filetag}.csv"
    local matrix_csv="${OUTDIR}/variant_line_matrix_${filetag}.csv"
    local tmpdir="${OUTDIR}/tmp_${filetag}"

    mkdir -p "$tmpdir"

    echo ""
    echo "============================================================"
    echo "[INFO] Processing: ${tag}"
    echo "============================================================"

    # -------------------------------------------------------------------------
    # STEP 1: Identify ancestor sample indices in the VCF
    # -------------------------------------------------------------------------

    echo "[INFO] Reading sample list from VCF..."

    mapfile -t all_samples < <(bcftools query -l "$vcf")
    n_samples=${#all_samples[@]}

    anc1_idx=-1
    anc2_idx=-1
    for i in "${!all_samples[@]}"; do
        s="${all_samples[$i]}"
        if [[ "$s" == "${ANC1}"* ]]; then anc1_idx=$i; fi
        if [[ "$s" == "${ANC2}"* ]]; then anc2_idx=$i; fi
    done

    if [[ $anc1_idx -lt 0 || $anc2_idx -lt 0 ]]; then
        echo "ERROR: Could not find both ancestors in VCF sample list." >&2
        echo "       ANC1 ($ANC1) index: $anc1_idx" >&2
        echo "       ANC2 ($ANC2) index: $anc2_idx" >&2
        echo "       Samples in VCF: ${all_samples[*]}" >&2
        exit 1
    fi

    echo "[INFO] Found ${ANC1} at index ${anc1_idx}"
    echo "[INFO] Found ${ANC2} at index ${anc2_idx}"
    echo "[INFO] Total samples: ${n_samples}"

    # -------------------------------------------------------------------------
    # STEP 2: Filter to sites where BOTH ancestors are hom-ref with GQ >= minGQ
    # -------------------------------------------------------------------------

    local anc_filtered="${tmpdir}/anc_filtered.vcf.gz"

    local anc_expr="( (GT[${anc1_idx}]=\"0/0\" || GT[${anc1_idx}]=\"0|0\") && GQ[${anc1_idx}]>=${ANC_MIN_GQ} ) && ( (GT[${anc2_idx}]=\"0/0\" || GT[${anc2_idx}]=\"0|0\") && GQ[${anc2_idx}]>=${ANC_MIN_GQ} )"

    echo "[INFO] Applying ancestor filter..."
    bcftools view -i "$anc_expr" "$vcf" -Oz -o "$anc_filtered"
    bcftools index -t "$anc_filtered"

    local n_before n_after
    n_before=$(bcftools view -H "$vcf"          | wc -l)
    n_after=$( bcftools view -H "$anc_filtered" | wc -l)
    echo "[INFO] Sites before ancestor filter : ${n_before}"
    echo "[INFO] Sites after ancestor filter  : ${n_after}"

    # -------------------------------------------------------------------------
    # STEP 3: Extract long-format table of hom-alt calls (non-ancestor only)
    # -------------------------------------------------------------------------
    # Output columns: CHROM  POS  REF  ALT  sample_full  line_key
    # line_key strips from first dot: HK_204.hifi_reads → HK_204

    echo "[INFO] Extracting hom-alt calls per sample..."

    local long_table="${tmpdir}/long_homalt.tsv.gz"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' \
        "$anc_filtered" \
    | awk -v OFS="\t" \
          -v anc1="${all_samples[$anc1_idx]}" \
          -v anc2="${all_samples[$anc2_idx]}" '
        function is_homalt(gt) {
            return (gt == "1/1" || gt == "1|1")
        }
        function short_key(samp,    tmp) {
            tmp = samp
            sub(/\..*/, "", tmp)
            return tmp
        }
        {
            chrom=$1; pos=$2; ref=$3; alt=$4
            for (i=5; i<=NF; i++) {
                split($i, a, "=")
                samp=a[1]; gt=a[2]
                if (samp == anc1 || samp == anc2) continue
                if (is_homalt(gt)) {
                    key = short_key(samp)
                    print chrom, pos, ref, alt, samp, key
                }
            }
        }
    ' | gzip -c > "$long_table"

    local n_calls
    n_calls=$(zcat "$long_table" | wc -l)
    echo "[INFO] Total hom-alt calls in non-ancestor samples: ${n_calls}"

    # -------------------------------------------------------------------------
    # STEP 4: Collapse to per-site summary
    # -------------------------------------------------------------------------
    # Output columns: CHROM  POS  REF  ALT  n_lines  line_list

    local per_site="${tmpdir}/per_site.tsv.gz"

    zcat "$long_table" \
    | sort -k1,1 -k2,2n -k3,3 -k4,4 -k6,6 \
    | awk -v OFS="\t" -v sep="$SEP" '
        BEGIN { prev=""; n=0; llist="" }

        function flush() {
            if (prev == "") return
            print prev_chrom, prev_pos, prev_ref, prev_alt, n, llist
        }

        {
            site = $1"\t"$2"\t"$3"\t"$4
            key  = $6

            if (site != prev && prev != "") {
                flush()
                delete seen
                n=0; llist=""
            }

            if (site != prev) {
                prev=site; prev_chrom=$1; prev_pos=$2; prev_ref=$3; prev_alt=$4
            }

            if (!(key in seen)) {
                seen[key]=1; n++
                llist = (llist == "") ? key : llist sep key
            }
        }

        END { flush() }
    ' | gzip -c > "$per_site"

    # -------------------------------------------------------------------------
    # STEP 5: Write contaminant_sites CSV
    # -------------------------------------------------------------------------
    # Sites where n_lines >= 2, sorted by n_lines descending

    {
        echo "CHROM,POS,REF,ALT,n_lines,line_list"
        zcat "$per_site" \
        | awk -F"\t" -v OFS="," '$5 >= 2 { print $1,$2,$3,$4,$5,$6 }' \
        | sort -t',' -k5,5rn
    } > "$contam_csv"

    local n_contam
    n_contam=$(tail -n +2 "$contam_csv" | wc -l)
    echo "[INFO] Contaminant/shared sites (n_lines >= 2): ${n_contam}"
    echo "[INFO] Written: ${contam_csv}"

    # -------------------------------------------------------------------------
    # STEP 6: Build the variant-sharing matrix
    # -------------------------------------------------------------------------

    echo "[INFO] Building variant-sharing matrix..."

    local lines_file="${tmpdir}/lines.txt"
    zcat "$per_site" \
    | awk -F"\t" -v sep="$SEP" '{
        n=split($6, a, sep)
        for(i=1;i<=n;i++) print a[i]
    }' \
    | sort -u > "$lines_file"

    local n_lines
    n_lines=$(wc -l < "$lines_file")
    echo "[INFO] Distinct non-ancestor lines: ${n_lines}"

    local per_site_txt="${tmpdir}/per_site.tsv"
    zcat "$per_site" > "$per_site_txt"

    awk -v FS="\t" -v OFS="," -v sep="$SEP" '
        FNR==NR {
            lines[++L] = $1
            next
        }
        {
            n = $5 + 0
            if (n <= 0) next
            split($6, a, sep)

            if (n == 1) {
                M[a[1], a[1]]++
            } else {
                for (i = 1; i <= n; i++) {
                    for (j = 1; j <= n; j++) {
                        if (i != j) M[a[i], a[j]]++
                    }
                }
            }
        }

        END {
            printf "Line"
            for (j = 1; j <= L; j++) printf "%s%s", OFS, lines[j]
            printf "\n"

            for (i = 1; i <= L; i++) {
                printf "%s", lines[i]
                for (j = 1; j <= L; j++) {
                    printf "%s%d", OFS, M[lines[i], lines[j]] + 0
                }
                printf "\n"
            }
        }
    ' "$lines_file" "$per_site_txt" > "$matrix_csv"

    echo "[INFO] Written: ${matrix_csv}"

    # -------------------------------------------------------------------------
    # STEP 7: Cleanup temporary files
    # -------------------------------------------------------------------------
    rm -rf "$tmpdir"

    echo "[INFO] Done: ${tag}"
}

# -----------------------------------------------------------------------------
# SECTION 5: Run for all four type × coverage combinations
# -----------------------------------------------------------------------------

for tag in snps.DP3 snps.DP10 indels.DP3 indels.DP10; do
    process_one "$tag"
done

# -----------------------------------------------------------------------------
# SECTION 6: Final summary
# -----------------------------------------------------------------------------

echo ""
echo "================================================================"
echo "[DONE] All done. Results in: ${OUTDIR}/"
echo "================================================================"
echo ""
echo "Contaminant site counts per dataset:"
for tag in snps.DP3 snps.DP10 indels.DP3 indels.DP10; do
    filetag="${tag//./_}"
    f="${OUTDIR}/contaminant_sites_${filetag}.csv"
    if [[ -f "$f" ]]; then
        n=$(tail -n +2 "$f" | wc -l)
        printf "  %-25s : %d shared sites\n" "$filetag" "$n"
    fi
done

echo ""
echo "Output files:"
for tag in snps.DP3 snps.DP10 indels.DP3 indels.DP10; do
    filetag="${tag//./_}"
    echo "  ${OUTDIR}/contaminant_sites_${filetag}.csv"
    echo "  ${OUTDIR}/variant_line_matrix_${filetag}.csv"
done
