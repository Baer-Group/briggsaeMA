#!/bin/sh

#SBATCH --job-name=variant_analysis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=120:00:00
#SBATCH --output=variant_analysis_%a.out
#SBATCH --account=juannanzhou --qos=juannanzhou
#SBATCH --array=0-3

cd /blue/baer/m.rifat/briggsae_MA_v2

#### Modules ####
module purge
module load bcftools

# Map array task ID to type/coverage combination
COMBOS=("snp_3x" "snp_10x" "indel_3x" "indel_10x")
TAG=${COMBOS[$SLURM_ARRAY_TASK_ID]}

TYPE=${TAG%_*}    # snp or indel
COV=${TAG#*_}     # 3x or 10x

INPUT_DIR="./filtered"
OUTPUT_DIR="./analysis_output"
LIST_SEP=";"

mkdir -p "${OUTPUT_DIR}"

b1="${INPUT_DIR}/B1.bp.GVCFs.${TYPE}.${COV}.filtered.biallelic.af_filtered.vcf"
b2="${INPUT_DIR}/B2.bp.GVCFs.${TYPE}.${COV}.filtered.biallelic.af_filtered.vcf"
tmpdir="${OUTPUT_DIR}/tmp_${TAG}"

mkdir -p "$tmpdir"

echo "========================================"
echo "Task ${SLURM_ARRAY_TASK_ID}: ${TAG}"
echo "B1: $(basename $b1)"
echo "B2: $(basename $b2)"
echo "========================================"

# =============================================================================
# FUNCTION: ancestor_filter
#   1. bgzip + index the input VCF
#   2. Identify all ANC* samples
#   3. Keep only sites where ALL ancestors are hom-ref (0/0 or 0|0) with GQ>=20
#   Outputs a bgzipped, indexed VCF into tmpdir
# =============================================================================
ancestor_filter() {
    local vcf="$1"
    local out_gz="$2"

    [[ -f "$vcf" ]] || { echo "ERROR: Missing input VCF: $vcf"; exit 1; }

    # bgzip + index
    local gz="${tmpdir}/$(basename ${vcf}).gz"
    echo "    Compressing: $(basename $vcf)"
    bcftools view -Oz -o "$gz" "$vcf"
    bcftools index -f "$gz"

    # Identify ANC sample indices
    mapfile -t all_samples < <(bcftools query -l "$gz")

    local anc_idx=()
    local anc_names=()
    for i in "${!all_samples[@]}"; do
        s="${all_samples[$i]}"
        if [[ "$s" == *ANC* ]]; then
            anc_idx+=("$i")
            anc_names+=("$s")
        fi
    done

    if (( ${#anc_idx[@]} == 0 )); then
        echo "    WARNING: No ANC samples found in $(basename $vcf) — skipping ancestor filter"
        cp -f "$gz" "$out_gz"
        bcftools index -f "$out_gz"
        rm -f "$gz" "${gz}.csi"
        return
    fi

    echo "    Ancestor samples found:"
    for j in "${!anc_idx[@]}"; do
        echo "      index ${anc_idx[$j]} -> ${anc_names[$j]}"
    done

    # Build bcftools filter expression (all ancestors must be hom-ref + GQ>=20)
    local expr=""
    for idx in "${anc_idx[@]}"; do
        expr+="( (GT[$idx]=\"0/0\" || GT[$idx]=\"0|0\") && GQ[$idx]>=20 ) && "
    done
    expr="${expr% && }"

    echo "    Filter expr: $expr"

    local before after
    before=$(bcftools view -H "$gz" | wc -l)

    bcftools view -i "$expr" -Oz -o "$out_gz" "$gz"
    bcftools index -f "$out_gz"

    after=$(bcftools view -H "$out_gz" | wc -l)
    echo "    Sites before ancestor filter : $before"
    echo "    Sites after ancestor filter  : $after"
    echo "    Sites removed                : $(( before - after ))"

    # Remove intermediate compressed copy
    rm -f "$gz" "${gz}.csi"
}

# =============================================================================
# FUNCTION: extract_long
#   Reads a bgzipped VCF, outputs gzipped long TSV of homalt variant calls.
#   Columns: CHROM  POS  REF  ALT  sample_full  line_key
#   line_key = numeric prefix before first underscore (e.g. 202)
# =============================================================================
extract_long() {
    local vcf_gz="$1"
    local out="$2"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' \
        "$vcf_gz" \
    | awk -v OFS="\t" '
        function is_homalt(gt) { return (gt == "1/1" || gt == "1|1") }
        function get_key(samp,   tmp) { tmp=samp; sub(/_.*/, "", tmp); return tmp }
        {
            chrom=$1; pos=$2; ref=$3; alt=$4
            for (i = 5; i <= NF; i++) {
                split($i, a, "=")
                samp=a[1]; gt=a[2]
                if (is_homalt(gt)) {
                    key=get_key(samp)
                    print chrom, pos, ref, alt, samp, key
                }
            }
        }
    ' | gzip -c > "$out"

    [[ -s "$out" ]] || { echo "ERROR: Empty extraction output: $out"; exit 1; }
}

# =============================================================================
# Step 0: Ancestor filter
# =============================================================================
echo ""
echo "[0/4] Applying ancestor filter..."

b1_anc="${tmpdir}/B1.anc_filtered.vcf.gz"
b2_anc="${tmpdir}/B2.anc_filtered.vcf.gz"

ancestor_filter "$b1" "$b1_anc"
ancestor_filter "$b2" "$b2_anc"

# =============================================================================
# Step 1: Extract long TSV tables from ancestor-filtered VCFs
# =============================================================================
echo ""
echo "[1/4] Extracting variant calls..."

extract_long "$b1_anc" "${tmpdir}/B1.long.tsv.gz"
extract_long "$b2_anc" "${tmpdir}/B2.long.tsv.gz"

# Ancestor-filtered VCFs no longer needed after extraction
rm -f "$b1_anc" "${b1_anc}.csi" "$b2_anc" "${b2_anc}.csi"

# =============================================================================
# Step 2: Combine, add site key, sort
# =============================================================================
echo ""
echo "[2/4] Combining and sorting..."

zcat "${tmpdir}/B1.long.tsv.gz" "${tmpdir}/B2.long.tsv.gz" \
| awk -v OFS="\t" '{ site=$1":"$2":"$3":"$4; print site,$1,$2,$3,$4,$5,$6 }' \
| sort -T "$tmpdir" --parallel=4 -k1,1 -k7,7 \
> "${tmpdir}/ALL.long.withsite.sorted.tsv"

gzip -c "${tmpdir}/ALL.long.withsite.sorted.tsv" \
> "${tmpdir}/ALL.long.withsite.sorted.tsv.gz"

# =============================================================================
# Step 3: Per-site collapse — builds per_site_lines and per_site_full in one pass
# =============================================================================
echo ""
echo "[3/4] Building per-site summaries..."

zcat "${tmpdir}/ALL.long.withsite.sorted.tsv.gz" \
| awk -v OFS="\t" -v sep="$LIST_SEP" \
    -v LINES_FILE="${tmpdir}/per_site_lines.tsv" \
    -v FULL_FILE="${tmpdir}/per_site_full.tsv" '
    BEGIN { prev_site=""; nkeys=0 }

    function flush() {
        if (prev_site == "") return
        print prev_site, prev_chrom, prev_pos, prev_ref, prev_alt, \
              nkeys, keys_list > LINES_FILE
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
        samples_list = (samples_list == "") ? samp : samples_list sep samp
        if (!(key in seen_key)) {
            seen_key[key]=1; nkeys++
            keys_list = (keys_list == "") ? key : keys_list sep key
        }
    }
    END { flush() }
'

gzip -c "${tmpdir}/per_site_lines.tsv" > "${tmpdir}/per_site_lines.tsv.gz"
gzip -c "${tmpdir}/per_site_full.tsv"  > "${tmpdir}/per_site.tsv.gz"

# =============================================================================
# Step 4a: Variant matrix
# =============================================================================
echo ""
echo "[4/4] Building variant matrix and contaminant outputs..."

# Build sorted line list
zcat "${tmpdir}/per_site_lines.tsv.gz" \
| awk -F"\t" -v sep="$LIST_SEP" '
    { n=split($7, a, sep); for (i=1; i<=n; i++) print a[i] }
' \
| sort -T "$tmpdir" -u -n \
> "${tmpdir}/lines.sorted"

# Square matrix
awk -v FS="\t" -v OFS="," -v sep="$LIST_SEP" '
    FNR == NR { lines[++L]=$1; next }
    {
        n=$6+0; if (n<=0) next
        split($7, a, sep)
        if (n==1) { M[a[1],a[1]]++ }
        else {
            for (i=1; i<=n; i++)
                for (j=i+1; j<=n; j++) { M[a[i],a[j]]++; M[a[j],a[i]]++ }
        }
    }
    END {
        printf "Line"
        for (j=1; j<=L; j++) printf "%s%s", OFS, lines[j]
        printf "\n"
        for (i=1; i<=L; i++) {
            printf "%s", lines[i]
            for (j=1; j<=L; j++) printf "%s%d", OFS, M[lines[i],lines[j]]+0
            printf "\n"
        }
    }
' "${tmpdir}/lines.sorted" "${tmpdir}/per_site_lines.tsv" \
> "${OUTPUT_DIR}/variant_matrix_${TAG}.csv"

# =============================================================================
# Step 4b: Contaminant outputs
# =============================================================================

# contaminant_sites.csv
# FIX: OFS must be passed via -v flag, NOT placed inline in the awk pattern
{
    echo "CHROM,POS,REF,ALT,n_keys,keys_list,samples_list"
    zcat "${tmpdir}/per_site.tsv.gz" \
    | awk -F"\t" -v OFS="," '($6+0)>=2 { print $2,$3,$4,$5,$6,$7,$8 }'
} > "${OUTPUT_DIR}/contaminant_sites_${TAG}.csv"

# contaminant site list
zcat "${tmpdir}/per_site.tsv.gz" \
| awk -F"\t" '($6+0)>=2 { print $1 }' \
| sort -T "$tmpdir" -k1,1 \
> "${tmpdir}/contam_sites.list"

# contaminant_events.csv
join -t $'\t' -1 1 -2 1 \
    "${tmpdir}/contam_sites.list" \
    <(sort -T "$tmpdir" -k1,1 "${tmpdir}/ALL.long.withsite.sorted.tsv") \
> "${tmpdir}/contam_long.tsv"

{
    echo "CHROM,POS,REF,ALT,key,sample_full,lane"
    awk -F"\t" -v OFS="," '
        {
            chrom=$2; pos=$3; ref=$4; alt=$5; samp=$6; key=$7
            lane="NA"
            if (samp ~ /_L001$/) lane="L001"
            else if (samp ~ /_L002$/) lane="L002"
            print chrom,pos,ref,alt,key,samp,lane
        }
    ' "${tmpdir}/contam_long.tsv"
} > "${OUTPUT_DIR}/contaminant_events_${TAG}.csv"

# summary_per_key.csv
zcat "${tmpdir}/per_site.tsv.gz" \
| awk -F"\t" -v OFS="\t" -v sep="$LIST_SEP" \
    '($6+0)>=2 { print $1,$7 }' \
> "${tmpdir}/contam_site_keys.tsv"

awk -F"\t" -v OFS="\t" -v sep="$LIST_SEP" '
    { n=split($2, a, sep); for (i=1; i<=n; i++) print $1, a[i] }
' "${tmpdir}/contam_site_keys.tsv" \
> "${tmpdir}/contam_site_key_pairs.tsv"

awk -F"\t" '{ print $2"\t"$1 }' "${tmpdir}/contam_site_key_pairs.tsv" \
| sort -T "$tmpdir" -u \
| awk -F"\t" '{ c[$1]++ } END { for (k in c) print k, c[k] }' \
| sort -T "$tmpdir" -k1,1 \
> "${tmpdir}/sites_per_key.tsv"

awk -F"\t" -v OFS="\t" '
    { site=$1; key=$2; keys[site]=keys[site] OFS key }
    END {
        for (s in keys) {
            line=substr(keys[s],2); n=split(line, a, OFS)
            for (i=1;i<=n;i++) for (j=1;j<=n;j++) if (i!=j) print a[i],a[j]
        }
    }
' "${tmpdir}/contam_site_key_pairs.tsv" \
| sort -T "$tmpdir" -u \
| awk -F"\t" '{ p[$1]++ } END { for (k in p) print k, p[k] }' \
| sort -T "$tmpdir" -k1,1 \
> "${tmpdir}/partners_per_key.tsv"

{
    echo "key,n_contaminant_sites,n_partner_keys"
    join -t $'\t' -a1 -a2 -e 0 -o 0,1.2,2.2 \
        "${tmpdir}/sites_per_key.tsv" \
        "${tmpdir}/partners_per_key.tsv" \
    | awk -F"\t" -v OFS="," '{ print $1,$2,$3 }'
} > "${OUTPUT_DIR}/summary_per_key_${TAG}.csv"

# =============================================================================
# Cleanup
# =============================================================================
rm -rf "$tmpdir"

echo ""
echo "========================================"
echo "Task ${SLURM_ARRAY_TASK_ID} (${TAG}) complete."
echo "Outputs written to: ${OUTPUT_DIR}/"
echo "  variant_matrix_${TAG}.csv"
echo "  contaminant_sites_${TAG}.csv"
echo "  contaminant_events_${TAG}.csv"
echo "  summary_per_key_${TAG}.csv"
echo "========================================"
