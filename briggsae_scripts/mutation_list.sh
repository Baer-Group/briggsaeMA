#!/bin/bash

#SBATCH --job-name=mutation_list
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=24:00:00
#SBATCH --output=mutation_list_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou

cd /blue/baer/m.rifat/briggsae_MA_v2

#### Modules ####
module purge
module load bcftools

# ===========================================================================
# PARAMETER: change COV to switch between 3x / 5x / 7x / 10x etc.
# ===========================================================================
COV="10x"

# ---------------------------------------------------------------------------
# Input VCFs  (derived from COV — only thing you need to change above)
# ---------------------------------------------------------------------------
INPUT_DIR="/blue/baer/m.rifat/briggsae_MA_v2/mutations"
B1_SNP="${INPUT_DIR}/B1.bp.GVCFs.snp.${COV}.mutation.vcf.gz"
B2_SNP="${INPUT_DIR}/B2.bp.GVCFs.snp.${COV}.mutation.vcf.gz"
B1_INDEL="${INPUT_DIR}/B1.bp.GVCFs.indel.${COV}.mutation.vcf.gz"
B2_INDEL="${INPUT_DIR}/B2.bp.GVCFs.indel.${COV}.mutation.vcf.gz"

OUTPUT_DIR="./analysis_output"
TMPDIR="${OUTPUT_DIR}/tmp_mutlist_${COV}"
FINAL_OUT="${OUTPUT_DIR}/mutation_list_${COV}.tsv"

mkdir -p "${OUTPUT_DIR}" "${TMPDIR}"

echo "========================================"
echo "mutation_list.sh  |  COV = ${COV}"
echo "========================================"

# ===========================================================================
# FUNCTION: extract_gt
#   Uses bcftools query to pull per-sample genotypes from a VCF.
#   Outputs a sorted, gzipped TSV:
#       CHROM  POS  REF  ALT  KEY  SHORT  GT
#
#   KEY   = numeric prefix for regular lines (e.g. 202),
#           or everything before _POOLRET for ANC samples (e.g. HK104_ANC_1)
#   SHORT = KEY_LANE  (e.g. 202_L001 or HK104_ANC_1_L001)
# ===========================================================================
extract_gt() {
    local vcf="$1"
    local out_gz="$2"   # final sorted .tsv.gz

    [[ -f "$vcf" ]] || { echo "ERROR: Missing VCF: $vcf"; exit 1; }

    # Files are already bgzipped — just index if needed, then query directly
    echo "  Indexing: $(basename $vcf)"
    bcftools index -f "$vcf"

    echo "  Extracting genotypes..."
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' "$vcf" \
    | awk -v OFS="\t" '
        function get_key(samp) {
            # Regular line: starts with digit -> numeric prefix only
            if (samp ~ /^[0-9]/) {
                key = samp; sub(/_.*/, "", key); return key
            }
            # ANC sample: take everything before _POOLRET
            key = samp
            if (match(key, /_POOLRET/)) key = substr(key, 1, RSTART-1)
            return key
        }
        function get_lane(samp) {
            if (match(samp, /_L[0-9]+$/)) return substr(samp, RSTART+1)
            return ""
        }
        {
            chrom=$1; pos=$2; ref=$3; alt=$4
            for (i = 5; i <= NF; i++) {
                split($i, a, "=")
                samp = a[1]; gt = a[2]
                key   = get_key(samp)
                lane  = get_lane(samp)
                short = (lane != "") ? key"_"lane : key
                print chrom, pos, ref, alt, key, short, gt
            }
        }
    ' \
    | sort -T "${TMPDIR}" --parallel=4 -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 \
    | gzip -c > "${out_gz}"

    echo "  Done: $(zcat ${out_gz} | wc -l) sample×site records"
}

# ===========================================================================
# FUNCTION: merge_and_filter
#   Full outer join on (CHROM, POS, REF, ALT, KEY) across B1 and B2.
#   Every site×sample seen in either replicate produces one output row,
#   as long as at least one replicate is hom-alt (1/1).
#
#   Column values per replicate:
#       1/1        → sample_L00X  (mutation confirmed)
#       0/0        → 0/0          (site present, ref call)
#       ./.        → missing      (site present, no call)
#       not in VCF → absent
#
#   Output columns (TSV, no header):
#       CHROM  POS  REF  ALT  KEY  B1  B2  TYPE
# ===========================================================================
merge_and_filter() {
    local b1_gz="$1"
    local b2_gz="$2"
    local type="$3"
    local out_tsv="$4"

    python3 - "$b1_gz" "$b2_gz" "$type" "$out_tsv" << 'PYEOF'
import sys, gzip

b1_path, b2_path, vtype, out_path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def is_homalt(gt): return gt in ("1/1", "1|1")
def is_missing(gt): return gt in ("./.", ".|.")

def fmt(gt, short):
    """Translate a raw genotype into the output column value."""
    if is_homalt(gt):  return short      # mutation present → sample name
    if is_missing(gt): return "missing"  # called but no genotype
    if gt in ("0/0", "0|0"): return "0/0"  # normalise phased ref to 0/0
    return gt                            # any other call

def is_mut_col(col):
    """True if the column value represents a confirmed mutation (sample name)."""
    return col not in ("0/0", "absent", "missing")

def records(path):
    """Yield (CHROM, POS_int, REF, ALT, KEY, SHORT, GT) from sorted gzipped TSV."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            yield p[0], int(p[1]), p[2], p[3], p[4], p[5], p[6]

def next_or_none(it):
    try: return next(it)
    except StopIteration: return None

# -----------------------------------------------------------------------
# Full outer join: walk both sorted iterators in lock-step.
# When only one side has a key, the other gets "absent".
# Only emit rows where at least one replicate is hom-alt (a mutation).
# -----------------------------------------------------------------------
written = 0
with open(out_path, "w") as out:
    iter1 = records(b1_path)
    iter2 = records(b2_path)

    r1 = next_or_none(iter1)
    r2 = next_or_none(iter2)

    while r1 is not None or r2 is not None:

        k1 = (r1[0], r1[1], r1[2], r1[3], r1[4]) if r1 is not None else None
        k2 = (r2[0], r2[1], r2[2], r2[3], r2[4]) if r2 is not None else None

        if k2 is None or (k1 is not None and k1 < k2):
            # Site×sample only in B1 → B2 is absent
            chrom, pos, ref, alt, key = k1
            b1_col = fmt(r1[6], r1[5])
            b2_col = "absent"
            r1 = next_or_none(iter1)

        elif k1 is None or k2 < k1:
            # Site×sample only in B2 → B1 is absent
            chrom, pos, ref, alt, key = k2
            b1_col = "absent"
            b2_col = fmt(r2[6], r2[5])
            r2 = next_or_none(iter2)

        else:
            # Both replicates have a call for this site×sample
            chrom, pos, ref, alt, key = k1
            b1_col = fmt(r1[6], r1[5])
            b2_col = fmt(r2[6], r2[5])
            r1 = next_or_none(iter1)
            r2 = next_or_none(iter2)

        # Only emit if at least one replicate shows the mutation
        if is_mut_col(b1_col) or is_mut_col(b2_col):
            out.write("\t".join([chrom, str(pos), ref, alt, key,
                                 b1_col, b2_col, vtype]) + "\n")
            written += 1

print(f"  [{vtype}] Mutations written: {written}", flush=True)
PYEOF
}


# ===========================================================================
# Main pipeline
# ===========================================================================

echo ""
echo "[1/3] Extracting genotype tables..."

extract_gt "$B1_SNP"   "${TMPDIR}/B1.snp.gt.tsv.gz"
extract_gt "$B2_SNP"   "${TMPDIR}/B2.snp.gt.tsv.gz"
extract_gt "$B1_INDEL" "${TMPDIR}/B1.indel.gt.tsv.gz"
extract_gt "$B2_INDEL" "${TMPDIR}/B2.indel.gt.tsv.gz"

echo ""
echo "[2/3] Merging replicates and applying filters..."

merge_and_filter \
    "${TMPDIR}/B1.snp.gt.tsv.gz" \
    "${TMPDIR}/B2.snp.gt.tsv.gz" \
    "snp" \
    "${TMPDIR}/mutations_snp.tsv"

merge_and_filter \
    "${TMPDIR}/B1.indel.gt.tsv.gz" \
    "${TMPDIR}/B2.indel.gt.tsv.gz" \
    "indel" \
    "${TMPDIR}/mutations_indel.tsv"

echo ""
echo "[3/3] Combining SNP + indel -> final output..."

{
    printf "CHROM\tPOS\tREF\tALT\tSample\tB1\tB2\tTYPE\n"
    cat "${TMPDIR}/mutations_snp.tsv" "${TMPDIR}/mutations_indel.tsv" \
    | sort --parallel=4 -T "${TMPDIR}" \
        -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 -k8,8
} > "${FINAL_OUT}"

TOTAL=$(( $(wc -l < "${FINAL_OUT}") - 1 ))
echo "  Total mutations: ${TOTAL}"

# ===========================================================================
# Cleanup
# ===========================================================================
rm -rf "${TMPDIR}"

echo ""
echo "========================================"
echo "Done. Output: ${FINAL_OUT}"
echo "========================================"
