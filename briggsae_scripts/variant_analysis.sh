#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Combined Variant Matrix + Contaminant Filter Script
#
# Input directory : ./filtered/
# Input filename  : B1/B2.bp.GVCFs.{snp/indel}.{3x/10x}.filtered.biallelic.af_filtered.vcf
#
# Output directory: ./analysis_output/
# Output files per combination (e.g. snp_3x):
#   variant_matrix_snp_3x.csv
#   contaminant_sites_snp_3x.csv
#   contaminant_events_snp_3x.csv
#   summary_per_key_snp_3x.csv
#
# Line key: numeric prefix before first underscore
#   e.g. 202_POOLRET91_S9_L001 -> 202
#
# Filters applied (homalt only — no hets in input):
#   Variant matrix : counts sites per line and shared sites between lines
#   Contaminant    : sites present as 1/1 in >=2 distinct biological keys
# ============================================================

module purge  >/dev/null 2>&1 || true
module load bcftools >/dev/null 2>&1 || true

# ---- user-tunable -----------------------------------------------------------
INPUT_DIR="./filtered"
OUTPUT_DIR="./analysis_output"
MODE="homalt"       # homalt only (no hets in input after previous filtering)
LIST_SEP=";"        # delimiter for lists inside CSV cells
# -----------------------------------------------------------------------------

die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing executable: $1"; }

need bcftools
need awk
need sort
need gzip
need zcat
need join

mkdir -p "${OUTPUT_DIR}"

# =============================================================================
# FUNCTION: extract_long
#   Reads an uncompressed VCF, outputs a gzipped long TSV of variant calls.
#   Columns: CHROM  POS  REF  ALT  sample_full  line_key
#
#   line_key  = numeric prefix before first underscore (e.g. 202)
#   key_style = "line"  -> strip everything after first underscore (matrix)
#             = "bio"   -> same rule, used for contaminant detection
# =============================================================================
extract_long() {
    local vcf="$1"
    local out="$2"

    [[ -f "$vcf" ]] || die "Missing input VCF: $vcf"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' \
        "$vcf" \
    | awk -v OFS="\t" '
        function is_homalt(gt) {
            return (gt == "1/1" || gt == "1|1")
        }
        function get_key(samp,   tmp) {
            tmp = samp
            sub(/_.*/, "", tmp)
            return tmp
        }
        {
            chrom=$1; pos=$2; ref=$3; alt=$4
            for (i = 5; i <= NF; i++) {
                split($i, a, "=")
                samp = a[1]; gt = a[2]
                if (is_homalt(gt)) {
                    key = get_key(samp)
                    print chrom, pos, ref, alt, samp, key
                }
            }
        }
    ' | gzip -c > "$out"

    [[ -s "$out" ]] || die "Extraction produced empty file: $out — check input VCF: $vcf"
}

# =============================================================================
# FUNCTION: build_variant_matrix
#   Reads per_site_lines.tsv.gz (already built) and writes the square matrix.
#   Diagonal  = sites unique to that line
#   Off-diag  = sites shared between two lines
# =============================================================================
build_variant_matrix() {
    local per_site_gz="$1"   # input: per_site_lines.tsv.gz
    local tmpdir="$2"
    local out_csv="$3"

    # Build sorted line list
    zcat "$per_site_gz" \
    | awk -F"\t" -v sep="$LIST_SEP" '
        { n = split($7, a, sep); for (i = 1; i <= n; i++) print a[i] }
    ' \
    | sort -T "$tmpdir" -u -n \
    > "${tmpdir}/lines.sorted"

    # Decompress per-site file for awk two-file pass
    zcat "$per_site_gz" > "${tmpdir}/per_site_lines.tsv"

    awk -v FS="\t" -v OFS="," -v sep="$LIST_SEP" '
        FNR == NR {
            lines[++L] = $1
            next
        }
        {
            n = $6 + 0
            if (n <= 0) next
            split($7, a, sep)
            if (n == 1) {
                M[a[1], a[1]]++
            } else {
                for (i = 1; i <= n; i++) {
                    for (j = i + 1; j <= n; j++) {
                        M[a[i], a[j]]++
                        M[a[j], a[i]]++
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
                    val = M[lines[i], lines[j]] + 0
                    printf "%s%d", OFS, val
                }
                printf "\n"
            }
        }
    ' "${tmpdir}/lines.sorted" "${tmpdir}/per_site_lines.tsv" > "$out_csv"
}

# =============================================================================
# FUNCTION: build_contaminant_outputs
#   Reads per_site.tsv.gz (already built) and produces:
#     contaminant_sites.csv
#     contaminant_events.csv
#     summary_per_key.csv
# =============================================================================
build_contaminant_outputs() {
    local per_site_gz="$1"   # input: per_site.tsv.gz (has keys + sample lists)
    local long_sorted_tsv="$2"  # input: ALL.long.withsite.sorted.tsv (uncompressed)
    local tmpdir="$3"
    local out_sites="$4"
    local out_events="$5"
    local out_summary="$6"

    # --- contaminant_sites.csv ---
    {
        echo "CHROM,POS,REF,ALT,n_keys,keys_list,samples_list"
        zcat "$per_site_gz" \
        | awk -F"\t" 'OFS="," ($6+0) >= 2 { print $2,$3,$4,$5,$6,$7,$8 }'
    } > "$out_sites"

    # --- contaminant site list (for joining) ---
    zcat "$per_site_gz" \
    | awk -F"\t" '($6+0) >= 2 { print $1 }' \
    | sort -T "$tmpdir" -k1,1 \
    > "${tmpdir}/contam_sites.list"

    # --- contaminant_events.csv ---
    join -t $'\t' -1 1 -2 1 \
        "${tmpdir}/contam_sites.list" \
        <(sort -T "$tmpdir" -k1,1 "$long_sorted_tsv") \
    > "${tmpdir}/contam_long.tsv"

    {
        echo "CHROM,POS,REF,ALT,key,sample_full,lane"
        awk -F"\t" -v OFS="," '
            {
                chrom=$2; pos=$3; ref=$4; alt=$5; samp=$6; key=$7
                lane = "NA"
                if (samp ~ /_L001$/) lane = "L001"
                else if (samp ~ /_L002$/) lane = "L002"
                print chrom, pos, ref, alt, key, samp, lane
            }
        ' "${tmpdir}/contam_long.tsv"
    } > "$out_events"

    # --- summary_per_key.csv ---
    # Build site->keys pairs for contaminant sites only
    zcat "$per_site_gz" \
    | awk -F"\t" -v OFS="\t" -v sep="$LIST_SEP" \
        '($6+0) >= 2 { print $1, $7 }' \
    > "${tmpdir}/contam_site_keys.tsv"

    # Expand to (site, key) pairs
    awk -F"\t" -v OFS="\t" -v sep="$LIST_SEP" '
        {
            site = $1; keys = $2
            n = split(keys, a, sep)
            for (i = 1; i <= n; i++) print site, a[i]
        }
    ' "${tmpdir}/contam_site_keys.tsv" \
    > "${tmpdir}/contam_site_key_pairs.tsv"

    # Sites per key
    awk -F"\t" '{ print $2"\t"$1 }' "${tmpdir}/contam_site_key_pairs.tsv" \
    | sort -T "$tmpdir" -u \
    | awk -F"\t" '{ k=$1; c[k]++ } END { for (k in c) print k, c[k] }' \
    | sort -T "$tmpdir" -k1,1 \
    > "${tmpdir}/sites_per_key.tsv"

    # Partner keys per key
    awk -F"\t" -v OFS="\t" '
        {
            site=$1; key=$2
            keys[site] = keys[site] OFS key
        }
        END {
            for (s in keys) {
                line = substr(keys[s], 2)
                n = split(line, a, OFS)
                for (i = 1; i <= n; i++)
                    for (j = 1; j <= n; j++)
                        if (i != j) print a[i], a[j]
            }
        }
    ' "${tmpdir}/contam_site_key_pairs.tsv" \
    | sort -T "$tmpdir" -u \
    | awk -F"\t" '{ k=$1; p[k]++ } END { for (k in p) print k, p[k] }' \
    | sort -T "$tmpdir" -k1,1 \
    > "${tmpdir}/partners_per_key.tsv"

    {
        echo "key,n_contaminant_sites,n_partner_keys"
        join -t $'\t' -a1 -a2 -e 0 -o 0,1.2,2.2 \
            "${tmpdir}/sites_per_key.tsv" \
            "${tmpdir}/partners_per_key.tsv" \
        | awk -F"\t" -v OFS="," '{ print $1, $2, $3 }'
    } > "$out_summary"
}

# =============================================================================
# MAIN: loop over all type x coverage combinations
# =============================================================================
echo "========================================"
echo "Combined Variant Matrix + Contaminant Filter"
echo "Input  : ${INPUT_DIR}"
echo "Output : ${OUTPUT_DIR}"
echo "Mode   : ${MODE} (homozygous ALT only)"
echo "========================================"

for type in snp indel; do
    for cov in 3x 10x; do

        tag="${type}_${cov}"
        b1="${INPUT_DIR}/B1.bp.GVCFs.${type}.${cov}.filtered.biallelic.af_filtered.vcf"
        b2="${INPUT_DIR}/B2.bp.GVCFs.${type}.${cov}.filtered.biallelic.af_filtered.vcf"
        tmpdir="${OUTPUT_DIR}/tmp_${tag}"

        mkdir -p "$tmpdir"

        echo ""
        echo "--> Processing: ${tag}"
        echo "    B1: $(basename $b1)"
        echo "    B2: $(basename $b2)"

        # ------------------------------------------------------------------
        # Step 1: Extract long TSV tables from each batch VCF
        # ------------------------------------------------------------------
        echo "    [1/4] Extracting variant calls..."
        extract_long "$b1" "${tmpdir}/B1.long.tsv.gz"
        extract_long "$b2" "${tmpdir}/B2.long.tsv.gz"

        # ------------------------------------------------------------------
        # Step 2: Combine both batches, add site key, sort
        # ------------------------------------------------------------------
        echo "    [2/4] Combining and sorting..."
        zcat "${tmpdir}/B1.long.tsv.gz" "${tmpdir}/B2.long.tsv.gz" \
        | awk -v OFS="\t" '{ site=$1":"$2":"$3":"$4; print site,$1,$2,$3,$4,$5,$6 }' \
        | sort -T "$tmpdir" -k1,1 -k7,7 \
        > "${tmpdir}/ALL.long.withsite.sorted.tsv"

        # Also gzip a copy for the per-site collapse step
        gzip -c "${tmpdir}/ALL.long.withsite.sorted.tsv" \
        > "${tmpdir}/ALL.long.withsite.sorted.tsv.gz"

        # ------------------------------------------------------------------
        # Step 3: Collapse to per-site unique line/key lists
        #   per_site_lines.tsv.gz : for variant matrix (col7 = line list)
        #   per_site.tsv.gz       : for contaminant (col7 = keys, col8 = samples)
        #   Both use the same key definition — built in one pass
        # ------------------------------------------------------------------
        echo "    [3/4] Building per-site summaries..."
        zcat "${tmpdir}/ALL.long.withsite.sorted.tsv.gz" \
        | awk -v OFS="\t" -v sep="$LIST_SEP" '
            BEGIN { prev_site=""; nkeys=0 }

            function flush() {
                if (prev_site == "") return
                # per_site_lines: site chrom pos ref alt n_lines lines_list
                print prev_site, prev_chrom, prev_pos, prev_ref, prev_alt, \
                      nkeys, keys_list > LINES_FILE
                # per_site: site chrom pos ref alt n_keys keys_list samples_list
                print prev_site, prev_chrom, prev_pos, prev_ref, prev_alt, \
                      nkeys, keys_list, samples_list > FULL_FILE
            }

            {
                site=$1; chrom=$2; pos=$3; ref=$4; alt=$5; samp=$6; key=$7

                if (site != prev_site && prev_site != "") {
                    flush()
                    delete seen_key
                    nkeys=0; keys_list=""; samples_list=""
                }
                if (site != prev_site) {
                    prev_site=site; prev_chrom=chrom; prev_pos=pos
                    prev_ref=ref; prev_alt=alt
                }

                # accumulate all sample names
                samples_list = (samples_list == "") ? samp : samples_list sep samp

                # unique key count
                if (!(key in seen_key)) {
                    seen_key[key] = 1
                    nkeys++
                    keys_list = (keys_list == "") ? key : keys_list sep key
                }
            }
            END { flush() }
        ' \
            LINES_FILE="${tmpdir}/per_site_lines.tsv" \
            FULL_FILE="${tmpdir}/per_site_full.tsv"

        gzip -c "${tmpdir}/per_site_lines.tsv" > "${tmpdir}/per_site_lines.tsv.gz"
        gzip -c "${tmpdir}/per_site_full.tsv"  > "${tmpdir}/per_site.tsv.gz"

        # ------------------------------------------------------------------
        # Step 4a: Variant matrix
        # ------------------------------------------------------------------
        echo "    [4/4] Building variant matrix and contaminant outputs..."
        build_variant_matrix \
            "${tmpdir}/per_site_lines.tsv.gz" \
            "$tmpdir" \
            "${OUTPUT_DIR}/variant_matrix_${tag}.csv"

        # ------------------------------------------------------------------
        # Step 4b: Contaminant outputs
        # ------------------------------------------------------------------
        build_contaminant_outputs \
            "${tmpdir}/per_site.tsv.gz" \
            "${tmpdir}/ALL.long.withsite.sorted.tsv" \
            "$tmpdir" \
            "${OUTPUT_DIR}/contaminant_sites_${tag}.csv" \
            "${OUTPUT_DIR}/contaminant_events_${tag}.csv" \
            "${OUTPUT_DIR}/summary_per_key_${tag}.csv"

        # ------------------------------------------------------------------
        # Cleanup temporary files
        # ------------------------------------------------------------------
        rm -rf "$tmpdir"

        echo "    Done: ${tag}"
        echo "      variant_matrix_${tag}.csv"
        echo "      contaminant_sites_${tag}.csv"
        echo "      contaminant_events_${tag}.csv"
        echo "      summary_per_key_${tag}.csv"

    done
done

echo ""
echo "========================================"
echo "All done. Outputs written to: ${OUTPUT_DIR}/"
echo "========================================"
