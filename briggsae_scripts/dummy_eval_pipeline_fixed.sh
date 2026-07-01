#!/usr/bin/env bash
#SBATCH --job-name=dummy_eval
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --time=120:00:00
#SBATCH --output=dummy_eval_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou

# =============================================================================
# C. briggsae Dummy Reference — Variant Recall Evaluation Pipeline
# =============================================================================
# PURPOSE:
#   Evaluate how well the DeepVariant pipeline recalls injected dummy variants.
#   For each injected site × each sample, classify the genotype as:
#       1/1  → recalled correctly
#       0/1  → heterozygous call  (failure to recall)
#       0/0  → homozygous ref     (false negative)
#       ./.  → missing/masked     (failure to recall)
#       absent → site not in VCF  (false negative, treated as 0/0)
#
# EVALUATION POINT:
#   Genotypes are assessed AFTER biallelic + DP masking but BEFORE the
#   AF filter, missing filter, and mutation filter. This is intentional:
#   - AF filter drops entire het sites → would erase 0/1 evidence
#   - Missing/mutation filters would hide FN and missing signals
#   - We want to see the raw per-sample genotype at each dummy site
#
# COORDINATE HANDLING:
#   SNP dummy : 1:1 coordinate mapping → evaluate directly
#   Indel dummy: cumulative offsets shift positions downstream of each indel
#               → ground truth VCF positions are lifted to dummy coords
#                 using indel_offset_table.csv before intersection
#
# OUTPUTS per run (SNP and indel, DP3 and DP10):
#   {tag}.biallelic.vcf.gz       — pipeline VCF used for evaluation
#   {tag}.gt_table.csv           — per-site × per-sample genotype table
#   {tag}.recall_summary.csv     — per-sample recall counts and rates
#   {tag}.by_length_summary.csv  — (indels only) recall by indel length bin
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# SECTION 1: User settings
# -----------------------------------------------------------------------------

# --- SNP dummy ---
SNP_WORKDIR="/orange/juannanzhou/Rifat_CB_raw/briggsae/dummy_snp"
SNP_REF="snp_dummy.fa"
SNP_GT_VCF="ground_truth_snps.fixed.vcf"   # in original = dummy coordinates

# --- Indel dummy ---
INDEL_WORKDIR="/orange/juannanzhou/Rifat_CB_raw/briggsae/dummy_indel"
INDEL_REF="indel_dummy.fa"
INDEL_GT_VCF="ground_truth_indels.fixed.vcf"    # in ORIGINAL coordinates
INDEL_OFFSET_TABLE="indel_offset_table.csv"  # produced by dummy_reference_briggsae.R

# --- Shared ---
THREADS=8
DP_LEVELS=(3 10)

# -----------------------------------------------------------------------------
# SECTION 2: Module setup
# -----------------------------------------------------------------------------

module purge || true
module load bcftools
module load samtools

# =============================================================================
# FUNCTION: run_dummy_pipeline
#   Runs GLnexus → BCF → VCF → biallelic for one dummy type
#   then calls the evaluation step
#
# Arguments:
#   $1 WORKDIR   — working directory containing gVCFs and reference
#   $2 REF       — dummy reference FASTA filename
#   $3 GT_VCF    — ground truth VCF (original coordinates)
#   $4 DUMMY_TYPE — "snp" or "indel"
#   $5 OFFSET_CSV — path to offset table (only used for indel, else "none")
# =============================================================================

run_dummy_pipeline() {

    local WORKDIR="$1"
    local REF="$2"
    local GT_VCF="$3"
    local DUMMY_TYPE="$4"
    local OFFSET_CSV="$5"
    local OUTDIR="${WORKDIR}/eval_out"

    echo ""
    echo "############################################################"
    echo "  DUMMY PIPELINE: ${DUMMY_TYPE^^}"
    echo "  Working dir: ${WORKDIR}"
    echo "############################################################"

    cd "$WORKDIR"
    mkdir -p "$OUTDIR"

    # Verify reference and ground truth exist
    [[ -f "$REF"    ]] || { echo "ERROR: Reference not found: $REF" >&2; exit 1; }
    [[ -f "$GT_VCF" ]] || { echo "ERROR: Ground truth VCF not found: $GT_VCF" >&2; exit 1; }

    # Index reference if needed
    [[ -f "${REF}.fai" ]] || samtools faidx "$REF"

    # -----------------------------------------------------------------
    # Step 1: Index gVCFs and build list
    # -----------------------------------------------------------------
    echo "[INFO] Indexing gVCFs..."
    GVCF_LIST="${OUTDIR}/gvcfs.list"
    ls *.dv.g.renamed.vcf.gz > "$GVCF_LIST"

    while read -r gvcf; do
        [[ -f "${gvcf}.tbi" ]] || bcftools index -t "$gvcf"
    done < "$GVCF_LIST"

    echo "[INFO] Found $(wc -l < "$GVCF_LIST") gVCF files."

    # -----------------------------------------------------------------
    # Step 2: GLnexus joint calling
    # -----------------------------------------------------------------
    module load glnexus

    BCF="${OUTDIR}/cohort.dummy_${DUMMY_TYPE}.bcf"

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

    # -----------------------------------------------------------------
    # Step 3: BCF → VCF
    # -----------------------------------------------------------------
    JOINT_VCF="${OUTDIR}/cohort.dummy_${DUMMY_TYPE}.vcf.gz"

    if [[ -f "$JOINT_VCF" && -f "${JOINT_VCF}.tbi" ]]; then
        echo "[INFO] Joint VCF already exists, skipping: $JOINT_VCF"
    else
        bcftools view "$BCF" -Oz -o "$JOINT_VCF"
        bcftools index -t "$JOINT_VCF"
        echo "[INFO] Joint VCF ready: $JOINT_VCF"
    fi

    # -----------------------------------------------------------------
    # Step 4: For each coverage depth → biallelic VCF → evaluate
    # -----------------------------------------------------------------
    for dp in "${DP_LEVELS[@]}"; do

        local tag="${DUMMY_TYPE}.DP${dp}"
        echo ""
        echo "------------------------------------------------------------"
        echo "[INFO] Processing: ${DUMMY_TYPE^^} | DP >= ${dp}"
        echo "------------------------------------------------------------"

        # Extract variant type, mask low-DP, keep biallelic
        # (same as main pipeline — normalization intentionally skipped)
        local BIALLELIC="${OUTDIR}/${tag}.biallelic.vcf.gz"

        bcftools view -v "${DUMMY_TYPE}s" "$JOINT_VCF" \
        | bcftools +setGT -- -t q -n . -i "FMT/DP<${dp}" \
        | bcftools view -m2 -M2 -Oz -o "$BIALLELIC"

        bcftools index -t "$BIALLELIC"

        local n_sites
        n_sites=$(bcftools view -H "$BIALLELIC" | wc -l)
        echo "[INFO] Sites in biallelic VCF: ${n_sites}"

        # -----------------------------------------------------------------
        # Step 5: Lift ground truth coordinates if indel
        # -----------------------------------------------------------------
        # For SNP dummy: original coords == dummy coords → use GT_VCF as-is
        # For indel dummy: lift original positions → dummy positions using
        # the cumulative offset table, then sort and bgzip for bcftools
        # -----------------------------------------------------------------

        local EVAL_GT_VCF
        if [[ "$DUMMY_TYPE" == "indel" ]]; then

            EVAL_GT_VCF="${OUTDIR}/ground_truth_indels.dummy_coords.vcf.gz"

            echo "[INFO] Lifting ground truth indel positions to dummy coordinates..."
            module load python/3.10

python3 << PYEOF
import csv, gzip

offset_csv = "${OFFSET_CSV}"
gt_vcf     = "${GT_VCF}"
out_vcf    = "${EVAL_GT_VCF%.gz}"   # write uncompressed first

# Build offset table: for each original position, what is dummy position?
# The offset table has rows: CHROM, ORIG_POS, DUMMY_POS, INDEL_TYPE, INDEL_LEN, CUM_OFFSET
# Strategy: for a given original position, find the largest ORIG_POS in the
# table that is <= query position, then apply that row's CUM_OFFSET.
from collections import defaultdict

offsets_by_chr = defaultdict(list)  # chr -> sorted list of (orig_pos, cum_offset)

with open(offset_csv) as f:
    reader = csv.DictReader(f)
    for row in reader:
        chrom    = row['CHROM']
        orig_pos = int(row['ORIG_POS'])
        cum_off  = int(row['CUM_OFFSET'])
        offsets_by_chr[chrom].append((orig_pos, cum_off))

# Sort each chromosome's list by orig_pos
for chrom in offsets_by_chr:
    offsets_by_chr[chrom].sort()

import bisect

def get_dummy_pos(chrom, orig_pos):
    """Return dummy coordinate for a given original coordinate."""
    entries = offsets_by_chr.get(chrom, [])
    if not entries:
        return orig_pos
    # Find the last indel at or before orig_pos
    positions = [e[0] for e in entries]
    idx = bisect.bisect_left(positions, orig_pos) - 1
    if idx < 0:
        return orig_pos   # no indels before this position, no offset yet
    cum_offset = entries[idx][1]
    return orig_pos + cum_offset

with open(gt_vcf) as fin, open(out_vcf, 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        fields = line.strip().split('\t')
        chrom    = fields[0]
        orig_pos = int(fields[1])
        dummy_pos = get_dummy_pos(chrom, orig_pos)
        fields[1] = str(dummy_pos)
        fout.write('\t'.join(fields) + '\n')

print(f"  Lifted ground truth to dummy coordinates: {out_vcf}")
PYEOF

            module unload python/3.10

            # Sort, bgzip, index the lifted VCF
            bcftools sort "${EVAL_GT_VCF%.gz}" -Oz -o "$EVAL_GT_VCF"
            bcftools index -t "$EVAL_GT_VCF"
            rm -f "${EVAL_GT_VCF%.gz}"

        else
            # SNP dummy — direct use, just bgzip + index if needed
            EVAL_GT_VCF="${OUTDIR}/ground_truth_snps.vcf.gz"
            if [[ ! -f "$EVAL_GT_VCF" ]]; then
                bcftools sort "$GT_VCF" -Oz -o "$EVAL_GT_VCF"
                bcftools index -t "$EVAL_GT_VCF"
            fi
        fi

        # -----------------------------------------------------------------
        # Step 6: Evaluate genotypes at ground truth positions
        # -----------------------------------------------------------------
        echo "[INFO] Evaluating genotypes at dummy sites..."
        module load python/3.10

        local GT_TABLE="${OUTDIR}/${tag}.gt_table.csv"
        local RECALL_SUMMARY="${OUTDIR}/${tag}.recall_summary.csv"
        local LEN_SUMMARY="${OUTDIR}/${tag}.by_length_summary.csv"

python3 << PYEOF
import gzip, csv
from collections import defaultdict

biallelic_vcf = "${BIALLELIC}"
gt_vcf        = "${EVAL_GT_VCF}"
gt_table_out  = "${GT_TABLE}"
recall_out    = "${RECALL_SUMMARY}"
len_out       = "${LEN_SUMMARY}"
dummy_type    = "${DUMMY_TYPE}"

# ------------------------------------------------------------------
# 1. Load ground truth sites → dict keyed by (CHROM, POS, REF, ALT)
#    Also store INFO field (has SIMTYPE, INDELSIZE)
# ------------------------------------------------------------------
gt_sites = {}   # (chrom, pos) -> (ref, alt, info_str)
with gzip.open(gt_vcf, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        key = (fields[0], int(fields[1]))
        gt_sites[key] = (fields[3], fields[4], fields[7])  # ref, alt, info

print(f"  Ground truth sites loaded: {len(gt_sites)}")

# ------------------------------------------------------------------
# 2. Load pipeline VCF → index by (CHROM, POS, REF, ALT) ->
#    dict of sample_name -> GT string
# ------------------------------------------------------------------
vcf_calls  = {}   # (chrom, pos, ref, alt) -> {sample: gt}
samples    = []

with gzip.open(biallelic_vcf, 'rt') as f:
    for line in f:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples = line.strip().split('\t')[9:]
            continue
        fields = line.strip().split('\t')
        chrom, pos = fields[0], int(fields[1])
        fmt    = fields[8].split(':')
        gt_idx = fmt.index('GT')
        key    = (chrom, pos)
        vcf_calls[key] = {}
        for i, samp in enumerate(samples):
            gt = fields[9+i].split(':')[gt_idx].replace('|','/')
            vcf_calls[key][samp] = gt

print(f"  Pipeline VCF sites loaded: {len(vcf_calls)}")
print(f"  Samples: {len(samples)}")

# ------------------------------------------------------------------
# 3. Classify per-site × per-sample
# ------------------------------------------------------------------
def classify_gt(gt):
    """Return category string for a genotype."""
    if gt in ('1/1','1|1'):       return '1/1'
    elif gt in ('0/1','0|1',
                '1/0','1|0'):     return '0/1'
    elif gt in ('0/0','0|0'):     return '0/0'
    elif '.' in gt:               return 'missing'
    else:                         return 'other'

def parse_indel_size(info_str):
    """Extract INDELSIZE from INFO field."""
    for f in info_str.split(';'):
        if f.startswith('INDELSIZE='):
            try: return int(f.split('=')[1])
            except: return None
    return None

def size_bin(size):
    if size is None:    return 'unknown'
    if size == 1:       return '1bp'
    if size == 2:       return '2bp'
    if 3  <= size <= 5: return '3-5bp'
    if 6  <= size <= 10:return '6-10bp'
    if 11 <= size <= 20:return '11-20bp'
    return '>20bp'

# Per-sample counters
sample_counts = {s: {'1/1':0,'0/1':0,'0/0':0,'missing':0,'other':0,'total':0}
                 for s in samples}

# Rows for GT table
gt_rows = []

# Count evaluable vs skipped sites
n_evaluable = 0
n_skipped   = 0

for (chrom, pos), (ref, alt, info) in gt_sites.items():
    indel_size = parse_indel_size(info) if dummy_type == 'indel' else None
    bin_label  = size_bin(indel_size)   if dummy_type == 'indel' else 'N/A'

    # Only evaluate sites present in the pipeline VCF.
    # Absent sites are unevaluable (spanned deletions, gVCF reference blocks,
    # or other dummy reference design artifacts) and are excluded from the
    # denominator to avoid inflating FN/FTR rates.
    if (chrom, pos) not in vcf_calls:
        n_skipped += 1
        continue

    n_evaluable += 1
    site_gts = vcf_calls[(chrom, pos)]

    for samp in samples:
        raw_gt = site_gts.get(samp, '0/0')
        cat    = classify_gt(raw_gt)
        sample_counts[samp][cat]   += 1
        sample_counts[samp]['total'] += 1
        gt_rows.append({
            'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
            'INFO': info, 'INDEL_BIN': bin_label,
            'SAMPLE': samp, 'GT': raw_gt, 'CATEGORY': cat
        })

print(f"  Evaluable sites : {n_evaluable} / {len(gt_sites)} total")
print(f"  Skipped (absent): {n_skipped} — unevaluable due to dummy reference design")

# ------------------------------------------------------------------
# 4. Write GT table
# ------------------------------------------------------------------
with open(gt_table_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=[
        'CHROM','POS','REF','ALT','INFO','INDEL_BIN',
        'SAMPLE','GT','CATEGORY'], delimiter=',')
    w.writeheader()
    w.writerows(gt_rows)
print(f"  GT table written: {gt_table_out}")

# ------------------------------------------------------------------
# 5. Write per-sample recall summary
# ------------------------------------------------------------------
with open(recall_out, 'w', newline='') as f:
    fieldnames = ['SAMPLE','Total','Recalled_1/1','Het_0/1','FN_0/0',
                  'Missing','Other','Recall_rate','FN_rate','FTR_rate']
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter=',')
    w.writeheader()
    for samp in samples:
        c   = sample_counts[samp]
        tot = c['total']
        if tot == 0: continue
        recalled = c['1/1']
        het      = c['0/1']
        fn       = c['0/0']
        miss     = c['missing']
        other    = c['other']
        ftr      = het + fn + miss + other   # failure to recall = anything != 1/1
        w.writerow({
            'SAMPLE':        samp,
            'Total':         tot,
            'Recalled_1/1':  recalled,
            'Het_0/1':       het,
            'FN_0/0':        fn,
            'Missing':       miss,
            'Other':         other,
            'Recall_rate':   f"{100*recalled/tot:.2f}%",
            'FN_rate':       f"{100*fn/tot:.2f}%",
            'FTR_rate':      f"{100*ftr/tot:.2f}%",
        })
print(f"  Recall summary written: {recall_out}")

# ------------------------------------------------------------------
# 6. Write indel-length-bin summary (indel only)
# ------------------------------------------------------------------
if dummy_type == 'indel':
    bin_counts = defaultdict(lambda: {'1/1':0,'0/1':0,'0/0':0,'missing':0,'total':0})
    for row in gt_rows:
        b = row['INDEL_BIN']
        bin_counts[b][row['CATEGORY']] += 1
        bin_counts[b]['total'] += 1

    bin_order = ['1bp','2bp','3-5bp','6-10bp','11-20bp','>20bp','unknown']
    with open(len_out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'INDEL_BIN','Total','Recalled_1/1','Het_0/1',
            'FN_0/0','Missing','Recall_rate','FN_rate','FTR_rate'],
            delimiter=',')
        w.writeheader()
        for b in bin_order:
            if b not in bin_counts: continue
            c   = bin_counts[b]
            tot = c['total']
            if tot == 0: continue
            recalled = c['1/1']
            ftr = tot - recalled
            w.writerow({
                'INDEL_BIN':    b,
                'Total':        tot,
                'Recalled_1/1': recalled,
                'Het_0/1':      c['0/1'],
                'FN_0/0':       c['0/0'],
                'Missing':      c['missing'],
                'Recall_rate':  f"{100*recalled/tot:.2f}%",
                'FN_rate':      f"{100*c['0/0']/tot:.2f}%",
                'FTR_rate':     f"{100*ftr/tot:.2f}%",
            })
    print(f"  Indel length summary written: {len_out}")

# ------------------------------------------------------------------
# 7. Print console summary
# ------------------------------------------------------------------
print(f"\n  === RECALL SUMMARY ===")
print(f"  {'SAMPLE':<40}  {'Total':>6}  {'1/1':>6}  {'0/1':>6}  {'0/0':>6}  "
      f"{'Missing':>8}  {'Recall%':>8}  {'FN%':>6}  {'FTR%':>6}")
print("  " + "-"*100)
all_tot = all_rec = all_het = all_fn = all_miss = 0
for samp in samples:
    c = sample_counts[samp]
    tot = c['total']
    if tot == 0: continue
    all_tot  += tot;      all_rec  += c['1/1']
    all_het  += c['0/1']; all_fn   += c['0/0']
    all_miss += c['missing']
    ftr = tot - c['1/1']
    print(f"  {samp:<40}  {tot:>6}  {c['1/1']:>6}  {c['0/1']:>6}  {c['0/0']:>6}  "
          f"{c['missing']:>8}  {100*c['1/1']/tot:>7.2f}%  "
          f"{100*c['0/0']/tot:>5.2f}%  {100*ftr/tot:>5.2f}%")

print("  " + "-"*100)
if all_tot > 0:
    all_ftr = all_tot - all_rec
    print(f"  {'OVERALL':<40}  {all_tot:>6}  {all_rec:>6}  {all_het:>6}  {all_fn:>6}  "
          f"{all_miss:>8}  {100*all_rec/all_tot:>7.2f}%  "
          f"{100*all_fn/all_tot:>5.2f}%  {100*all_ftr/all_tot:>5.2f}%")
PYEOF

        module unload python/3.10
        echo "[INFO] Evaluation complete for ${tag}"
        echo "[INFO] Outputs:"
        echo "         ${GT_TABLE}"
        echo "         ${RECALL_SUMMARY}"
        [[ "$DUMMY_TYPE" == "indel" ]] && echo "         ${LEN_SUMMARY}"

    done   # end DP loop

    echo ""
    echo "[DONE] ${DUMMY_TYPE^^} dummy evaluation complete."
    echo "       Results in: ${OUTDIR}"
}

# =============================================================================
# MAIN: Run SNP dummy then indel dummy
# =============================================================================

run_dummy_pipeline \
    "$SNP_WORKDIR" \
    "$SNP_REF" \
    "$SNP_GT_VCF" \
    "snp" \
    "none"

run_dummy_pipeline \
    "$INDEL_WORKDIR" \
    "$INDEL_REF" \
    "$INDEL_GT_VCF" \
    "indel" \
    "${INDEL_WORKDIR}/${INDEL_OFFSET_TABLE}"

# =============================================================================
# FINAL SUMMARY
# =============================================================================
echo ""
echo "================================================================"
echo "[DONE] All dummy evaluations complete."
echo "================================================================"
echo ""
echo "SNP results:   ${SNP_WORKDIR}/eval_out/"
echo "Indel results: ${INDEL_WORKDIR}/eval_out/"
echo ""
echo "Key output files per run (snp.DP3, snp.DP10, indel.DP3, indel.DP10):"
echo "  *.gt_table.csv         — full per-site × per-sample genotype table"
echo "  *.recall_summary.csv   — per-sample: 1/1, 0/1, 0/0, missing counts + rates"
echo "  *.by_length_summary.csv — (indels) recall broken down by indel length bin"
