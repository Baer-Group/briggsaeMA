#!/usr/bin/env bash
#SBATCH --job-name=b1_dummy_eval
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --output=b1_dummy_eval_%j.out
#SBATCH --account=baer --qos=baer

# =============================================================================
# C. elegans Short-Read Dummy Reference — Variant Recall Evaluation Pipeline
# =============================================================================
# PURPOSE:
#   Evaluate how well the GATK short-read pipeline recalls injected dummy
#   variants. For each injected site × each sample, classify the genotype as:
#       1/1  → recalled correctly
#       0/1  → heterozygous call  (failure to recall)
#       0/0  → homozygous ref     (false negative)
#       ./.  → missing/masked     (failure to recall)
#       absent → unevaluable (spanned deletion or ref block) — EXCLUDED
#
# IMPORTANT — DUMMY REFERENCE DESIGN:
#   The dummy reference has each injected ALT allele baked in as REF.
#   The pipeline therefore reports REF=original_ALT, ALT=original_REF.
#   Lookup uses (CHROM, POS) only — REF/ALT are NOT used as keys.
#   Absent sites are excluded from the denominator (not counted as FN).
#
# COORDINATE HANDLING:
#   SNP dummy : 1:1 coordinate mapping → evaluate directly
#   Indel dummy: cumulative offsets shift positions downstream of each indel
#               → ground truth VCF positions lifted to dummy coords using
#                 indel_offset_table.csv (bisect_left to avoid self-offset)
#
# KNOWN ISSUES FIXED:
#   1. Missing INFO headers (MUTTYPE, INDELSIZE) in ground truth VCFs
#      → fixed inline before bgzip/index
#   2. gzip.open() required for bgzipped ground truth VCFs
#   3. Position-only key lookup — REF/ALT swapped due to dummy reference
#   4. bisect_left (not bisect_right) for indel coordinate lifting
#   5. Absent sites excluded from denominator (unevaluable, not FN)
#   6. GATK bundles incomplete conda Python — unload gatk before python/3.10
#
# OUTPUTS per run (SNP and indel, DP3 and DP10):
#   {tag}.biallelic.vcf       — pipeline VCF used for evaluation
#   {tag}.gt_table.csv        — per-site × per-sample genotype table
#   {tag}.recall_summary.csv  — per-sample recall counts and rates
#   {tag}.by_length_summary.csv — (indels only) recall by indel length bin
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# SECTION 1: User settings
# -----------------------------------------------------------------------------

SNP_WORKDIR="/orange/baer/briggsae_dummy"
SNP_REF="snp_dummy.fa"
SNP_JOINT="B1_dummy_snp.bp.GVCFs.vcf"           # already exists
SNP_GT_VCF="ground_truth_snps.vcf"     # original coordinates

INDEL_WORKDIR="/orange/baer/briggsae_dummy_indel"
INDEL_REF="indel_dummy.fa"
INDEL_JOINT="B1_dummy_indel.bp.GVCFs.vcf"       # already exists
INDEL_GT_VCF="ground_truth_indels.vcf" # original coordinates
INDEL_OFFSET_TABLE="indel_offset_table.csv"    # produced by dummy_reference R script

DP_LEVELS=(3 10)
JAVA_MEM="-Xmx8g"

# -----------------------------------------------------------------------------
# SECTION 2: Module setup
# NOTE: gatk and python/3.10 conflict on HiPerGator — gatk is loaded only
# for DP filtering steps and unloaded before any python/3.10 call.
# -----------------------------------------------------------------------------

module purge || true
module load bcftools
module load samtools
module load gatk/

# =============================================================================
# FUNCTION: run_dummy_pipeline
#   Runs DP filtering → biallelic → evaluation for one dummy type
#
# Arguments:
#   $1 WORKDIR    — working directory
#   $2 REF        — dummy reference FASTA filename
#   $3 JOINT_VCF  — already-generated joint VCF (from GenotypeGVCFs)
#   $4 GT_VCF     — ground truth VCF (plain, original coordinates)
#   $5 DUMMY_TYPE — "snp" or "indel"
#   $6 OFFSET_CSV — path to offset table (indel only, else "none")
# =============================================================================

run_dummy_pipeline() {

    local WORKDIR="$1"
    local REF="$2"
    local JOINT_VCF="$3"
    local GT_VCF="$4"
    local DUMMY_TYPE="$5"
    local OFFSET_CSV="$6"
    local OUTDIR="${WORKDIR}/eval_out"

    echo ""
    echo "############################################################"
    echo "  DUMMY PIPELINE: ${DUMMY_TYPE^^}"
    echo "  Working dir: ${WORKDIR}"
    echo "############################################################"

    cd "$WORKDIR"
    mkdir -p "$OUTDIR"

    [[ -f "$REF"       ]] || { echo "ERROR: Reference not found: $REF" >&2; exit 1; }
    [[ -f "$JOINT_VCF" ]] || { echo "ERROR: Joint VCF not found: $JOINT_VCF" >&2; exit 1; }
    [[ -f "$GT_VCF"    ]] || { echo "ERROR: Ground truth VCF not found: $GT_VCF" >&2; exit 1; }

    [[ -f "${REF}.fai" ]] || samtools faidx "$REF"

    # -----------------------------------------------------------------
    # Step 1: Fix ground truth VCF headers and bgzip+index
    # -----------------------------------------------------------------
    # Ground truth VCFs use MUTTYPE (SNP) or INDELSIZE (indel) in records
    # but these INFO fields are missing from the header, causing bcftools
    # to error. Fix by inserting the missing ##INFO lines before bgzipping.
    # -----------------------------------------------------------------
    echo "[INFO] Fixing ground truth VCF headers..."
    module unload gatk

    local GT_VCF_GZ="${OUTDIR}/ground_truth_${DUMMY_TYPE}s.fixed.vcf.gz"

    if [[ ! -f "$GT_VCF_GZ" ]]; then
        if [[ "$DUMMY_TYPE" == "snp" ]]; then
            sed 's|##INFO=<ID=SIMTYPE,Number=1,Type=String,Description="Simulated variant type">|##INFO=<ID=SIMTYPE,Number=1,Type=String,Description="Simulated variant type">\n##INFO=<ID=MUTTYPE,Number=1,Type=String,Description="Mutation type (e.g. T>C)">|' \
                "$GT_VCF" \
            | bcftools sort -Oz -o "$GT_VCF_GZ"
        else
            sed 's|##INFO=<ID=SIMTYPE,Number=1,Type=String,Description="Simulated variant type">|##INFO=<ID=SIMTYPE,Number=1,Type=String,Description="Simulated variant type">\n##INFO=<ID=INDELSIZE,Number=1,Type=Integer,Description="Size of the indel in bp">|' \
                "$GT_VCF" \
            | bcftools sort -Oz -o "$GT_VCF_GZ"
        fi
        bcftools index -t "$GT_VCF_GZ"
        echo "[INFO] Fixed ground truth VCF: $GT_VCF_GZ"
    else
        echo "[INFO] Fixed ground truth VCF already exists, skipping."
    fi

    module load gatk/

    # -----------------------------------------------------------------
    # Step 2: For each coverage depth → DP filter → biallelic → evaluate
    # -----------------------------------------------------------------
    for dp in "${DP_LEVELS[@]}"; do

        local tag="${DUMMY_TYPE}.DP${dp}"
        echo ""
        echo "------------------------------------------------------------"
        echo "[INFO] Processing: ${DUMMY_TYPE^^} | DP >= ${dp}"
        echo "------------------------------------------------------------"

        # GATK DP filtering
        local FILTERED="${OUTDIR}/${tag}.vcf"
        local NOCALL="${OUTDIR}/${tag}.nocall.vcf"
        local CLEAN="${OUTDIR}/${tag}.filtered.vcf"
        local BIALLELIC="${OUTDIR}/${tag}.biallelic.vcf.gz"

        echo "[INFO] Applying DP${dp} filter..."

        gatk --java-options "$JAVA_MEM" VariantFiltration \
            -R "$REF" -V "$JOINT_VCF" \
            --genotype-filter-name "DP" \
            --genotype-filter-expression "DP < ${dp}.0" \
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

        # Biallelic filter
        bcftools view --max-alleles 2 "$CLEAN" -Oz -o "$BIALLELIC"
        bcftools index -t "$BIALLELIC"

        local n_sites
        n_sites=$(bcftools view -H "$BIALLELIC" | wc -l)
        echo "[INFO] Sites in biallelic VCF: ${n_sites}"

        # -----------------------------------------------------------------
        # Step 3: Lift ground truth coordinates if indel
        # NOTE: Unload gatk before loading python/3.10 to avoid the GATK
        # conda Python hijacking the module (causes ModuleNotFoundError)
        # -----------------------------------------------------------------

        local EVAL_GT_VCF
        if [[ "$DUMMY_TYPE" == "indel" ]]; then

            EVAL_GT_VCF="${OUTDIR}/ground_truth_indels.dummy_coords.vcf.gz"

            echo "[INFO] Lifting ground truth indel positions to dummy coordinates..."
            module unload gatk
            module load python/3.10

python3 << PYEOF
import csv, gzip
from collections import defaultdict
import bisect

offset_csv = "${OFFSET_CSV}"
gt_vcf_gz  = "${GT_VCF_GZ}"
out_vcf    = "${EVAL_GT_VCF%.gz}"

offsets_by_chr = defaultdict(list)

with open(offset_csv) as f:
    reader = csv.DictReader(f)
    for row in reader:
        chrom    = row['CHROM']
        orig_pos = int(row['ORIG_POS'])
        cum_off  = int(row['CUM_OFFSET'])
        offsets_by_chr[chrom].append((orig_pos, cum_off))

for chrom in offsets_by_chr:
    offsets_by_chr[chrom].sort()

def get_dummy_pos(chrom, orig_pos):
    entries = offsets_by_chr.get(chrom, [])
    if not entries:
        return orig_pos
    positions = [e[0] for e in entries]
    # bisect_left: only apply offsets from indels STRICTLY before this position
    idx = bisect.bisect_left(positions, orig_pos) - 1
    if idx < 0:
        return orig_pos
    return orig_pos + entries[idx][1]

with gzip.open(gt_vcf_gz, 'rt') as fin, open(out_vcf, 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        fields = line.strip().split('\t')
        chrom     = fields[0]
        orig_pos  = int(fields[1])
        dummy_pos = get_dummy_pos(chrom, orig_pos)
        fields[1] = str(dummy_pos)
        fout.write('\t'.join(fields) + '\n')

print(f"  Lifted ground truth to dummy coordinates: {out_vcf}")
PYEOF

            module unload python/3.10
            module load gatk/

            bcftools sort "${EVAL_GT_VCF%.gz}" -Oz -o "$EVAL_GT_VCF"
            bcftools index -t "$EVAL_GT_VCF"
            rm -f "${EVAL_GT_VCF%.gz}"

        else
            EVAL_GT_VCF="${OUTDIR}/ground_truth_snps.vcf.gz"
            if [[ ! -f "$EVAL_GT_VCF" ]]; then
                cp "$GT_VCF_GZ"  "$EVAL_GT_VCF"
                cp "${GT_VCF_GZ}.tbi" "${EVAL_GT_VCF}.tbi"
            fi
        fi

        # -----------------------------------------------------------------
        # Step 4: Evaluate genotypes at ground truth positions
        # Unload gatk before python/3.10 to avoid module conflict
        # -----------------------------------------------------------------
        echo "[INFO] Evaluating genotypes at dummy sites..."
        module unload gatk
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
# 1. Load ground truth sites → dict keyed by (CHROM, POS) only
#    REF/ALT stored for output but NOT used as lookup keys.
#    Dummy reference has injected ALT baked in as REF, so pipeline
#    reports swapped REF/ALT — position-only avoids this mismatch.
# ------------------------------------------------------------------
gt_sites = {}   # (chrom, pos) -> (ref, alt, info_str)
with gzip.open(gt_vcf, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        key = (fields[0], int(fields[1]))
        gt_sites[key] = (fields[3], fields[4], fields[7])

print(f"  Ground truth sites loaded: {len(gt_sites)}")

# ------------------------------------------------------------------
# 2. Load pipeline VCF → index by (CHROM, POS) only
# ------------------------------------------------------------------
vcf_calls = {}
samples   = []

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
# 3. Classify genotypes
# ------------------------------------------------------------------
def classify_gt(gt):
    if gt in ('1/1','1|1'):                     return '1/1'
    elif gt in ('0/1','0|1','1/0','1|0'):       return '0/1'
    elif gt in ('0/0','0|0'):                   return '0/0'
    elif '.' in gt:                             return 'missing'
    else:                                       return 'other'

def parse_indel_size(info_str):
    for f in info_str.split(';'):
        if f.startswith('INDELSIZE='):
            try: return int(f.split('=')[1])
            except: return None
    return None

def size_bin(size):
    if size is None:     return 'unknown'
    if size == 1:        return '1bp'
    if size == 2:        return '2bp'
    if 3  <= size <= 5:  return '3-5bp'
    if 6  <= size <= 10: return '6-10bp'
    if 11 <= size <= 20: return '11-20bp'
    return '>20bp'

sample_counts = {s: {'1/1':0,'0/1':0,'0/0':0,'missing':0,'other':0,'total':0}
                 for s in samples}
gt_rows    = []
n_evaluable = 0
n_skipped   = 0

for (chrom, pos), (ref, alt, info) in gt_sites.items():
    indel_size = parse_indel_size(info) if dummy_type == 'indel' else None
    bin_label  = size_bin(indel_size)   if dummy_type == 'indel' else 'N/A'

    # Skip sites absent from pipeline VCF (unevaluable — spanned deletions,
    # ref blocks, or other dummy reference design artifacts)
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
        ftr      = het + fn + miss + other
        w.writerow({
            'SAMPLE':       samp,
            'Total':        tot,
            'Recalled_1/1': recalled,
            'Het_0/1':      het,
            'FN_0/0':       fn,
            'Missing':      miss,
            'Other':        other,
            'Recall_rate':  f"{100*recalled/tot:.2f}%",
            'FN_rate':      f"{100*fn/tot:.2f}%",
            'FTR_rate':     f"{100*ftr/tot:.2f}%",
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
        module load gatk/
        echo "[INFO] Evaluation complete for ${tag}"

        # Cleanup intermediates for this tag
        rm -f "$FILTERED" "${FILTERED}.idx" \
              "$NOCALL"   "${NOCALL}.idx" \
              "$CLEAN"    "${CLEAN}.idx"

    done   # end DP loop

    echo ""
    echo "[DONE] ${DUMMY_TYPE^^} dummy evaluation complete."
    echo "       Results in: ${OUTDIR}"
}

# =============================================================================
# MAIN
# =============================================================================

run_dummy_pipeline \
    "$SNP_WORKDIR" \
    "$SNP_REF" \
    "$SNP_JOINT" \
    "$SNP_GT_VCF" \
    "snp" \
    "none"

run_dummy_pipeline \
    "$INDEL_WORKDIR" \
    "$INDEL_REF" \
    "$INDEL_JOINT" \
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
echo "  *.gt_table.csv          — full per-site × per-sample genotype table"
echo "  *.recall_summary.csv    — per-sample: 1/1, 0/1, 0/0, missing + rates"
echo "  *.by_length_summary.csv — (indels) recall broken down by indel length bin"
