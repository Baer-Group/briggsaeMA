#!/usr/bin/env bash
#SBATCH --job-name=anc_dv_joint
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --output=anc_dv_joint_%j.out
#SBATCH --account=baer --qos=baer-b

# =============================================================================
# Ancestor Joint Calling Pipeline — SR DeepVariant + GLnexus + Count Summary
# =============================================================================
# WORKDIR: /blue/baer/m.rifat/briggsae_ancestors/Joint_ancestor_call
#
# INPUT (all in WORKDIR):
#   SR BAMs : HK104_ANC_1_POOLRET91_S18_L001.markduplicates.bam  (+ L002)
#             PB800_ANC__1_POOLRET91_S12_L001.markduplicates.bam  (+ L002)
#   LR gVCFs: HK_G0_1.hifi_reads.dv.g.renamed.vcf.gz
#             PB_G0_2.hifi_reads.dv.g.renamed.vcf.gz
#   REF     : 20250626_c_briggsae_Feb2020.genome.fa
#
# STEPS:
#   1. Run DeepVariant (WGS model) on 4 SR BAMs → per-sample gVCFs
#   2. Build unified gVCF list (4 SR + 2 LR) → all DeepVariant format
#   3. GLnexus joint calling → BCF → VCF
#   4. DP filter (3x and 10x) + biallelic + SNP/indel split
#   5. Python count summary — strict filter (all 6 = 0/0 or 1/1)
#      Indel size filter: |len(ALT)-len(REF)| <= 30
#
# OUTPUT: Joint_ancestor_call/summary/combined_count_summary_strict.csv
# =============================================================================

set -euo pipefail

WORKDIR="/blue/baer/m.rifat/briggsae_ancestors/Joint_ancestor_call"
REF="20250626_c_briggsae_Feb2020.genome.fa"
OUTDIR="${WORKDIR}/summary"
THREADS=8

# Fix /tmp exhaustion — use SLURM-private scratch for DeepVariant temp files
export TMPDIR="/scratch/local/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

cd "$WORKDIR"
mkdir -p "$OUTDIR"

module purge || true
module load bcftools
module load samtools

# =============================================================================
# STEP 1: DeepVariant WGS on 4 SR BAMs
# Skip samples whose gVCF already exists and is valid
# =============================================================================

echo ""
echo "[STEP 1] Running DeepVariant (WGS) on SR BAMs..."
module load deepvariant/1.10.0

for BAM in *.markduplicates.bam; do
    base="${BAM%.markduplicates.bam}"
    OUT_VCF="${base}.dv.vcf.gz"
    OUT_GVCF="${base}.dv.g.vcf.gz"

    # Safe skip: both outputs must exist and be valid BGZF
    if [[ -f "$OUT_GVCF" ]] && bcftools view "$OUT_GVCF" > /dev/null 2>&1; then
        echo "[SKIP] $base — gVCF already exists and is valid."
        continue
    fi

    # Remove any truncated output before rerunning
    rm -f "$OUT_VCF" "${OUT_VCF}.tbi" "$OUT_GVCF" "${OUT_GVCF}.tbi"

    echo "[RUN] DeepVariant on ${base}..."

    deepvariant \
        --model_type=WGS \
        --ref="$REF" \
        --reads="$BAM" \
        --output_vcf="$OUT_VCF" \
        --output_gvcf="$OUT_GVCF" \
        --num_shards="$THREADS"

    echo "[DONE] ${base}"
done

module unload deepvariant/1.10.0
echo "[STEP 1] All SR gVCFs ready."

# =============================================================================
# STEP 2: Build unified gVCF list (4 SR DeepVariant + 2 LR DeepVariant)
# =============================================================================

echo ""
echo "[STEP 2] Building unified gVCF list..."

GVCF_LIST="${OUTDIR}/all_gvcfs.list"

# SR: newly generated DeepVariant gVCFs
ls HK104_ANC_1_POOLRET91_S18_L001.dv.g.vcf.gz \
   HK104_ANC_1_POOLRET91_S18_L002.dv.g.vcf.gz \
   PB800_ANC__1_POOLRET91_S12_L001.dv.g.vcf.gz \
   PB800_ANC__1_POOLRET91_S12_L002.dv.g.vcf.gz \
   HK_G0_1.hifi_reads.dv.g.renamed.vcf.gz \
   PB_G0_2.hifi_reads.dv.g.renamed.vcf.gz > "$GVCF_LIST"

echo "[INFO] gVCF list:"
cat "$GVCF_LIST"

# Index any unindexed gVCFs
while read -r gvcf; do
    if [[ ! -f "${gvcf}.tbi" ]]; then
        echo "[INFO] Indexing ${gvcf}..."
        bcftools index -t "$gvcf"
    fi
done < "$GVCF_LIST"

echo "[STEP 2] Done — $(wc -l < "$GVCF_LIST") gVCFs."

# =============================================================================
# STEP 3: GLnexus joint calling
# =============================================================================

echo ""
echo "[STEP 3] GLnexus joint calling..."

module load glnexus

BCF="${OUTDIR}/cohort.joint.bcf"
JOINT_VCF="${OUTDIR}/cohort.joint.vcf.gz"

if [[ -f "$JOINT_VCF" && -f "${JOINT_VCF}.tbi" ]]; then
    echo "[SKIP] Joint VCF already exists: $JOINT_VCF"
else
    rm -rf GLnexus.DB

    glnexus_cli --config DeepVariantWGS \
                --list "$GVCF_LIST" \
        | bcftools view -Ob -o "$BCF"

    bcftools index "$BCF"

    bcftools view "$BCF" -Oz -o "$JOINT_VCF"
    bcftools index -t "$JOINT_VCF"

    echo "[STEP 3] Joint VCF ready: $JOINT_VCF"

    # Quick sanity check
    echo "[INFO] Sample names in joint VCF:"
    bcftools query -l "$JOINT_VCF"

    echo "[INFO] Spot check — first 5 variant sites:"
    bcftools view -H "$JOINT_VCF" | head -5 | cut -f1-9 || true
fi

module unload glnexus

# =============================================================================
# STEP 4: DP filter + biallelic + SNP/indel split
# =============================================================================

echo ""
echo "[STEP 4] Applying DP filters and splitting by variant type..."

for dp in 3 5 7 10; do
    for vtype in snps indels; do
        tag="${vtype}.DP${dp}"
        out="${OUTDIR}/${tag}.biallelic.vcf.gz"

        if [[ -f "$out" && -f "${out}.tbi" ]]; then
            echo "[SKIP] Already exists: $out"
            continue
        fi

        echo "[INFO] Processing: ${vtype} | DP >= ${dp}..."

        # Use intermediate files to avoid SIGPIPE with set -euo pipefail
        tmp1="${OUTDIR}/${tag}.tmp1.vcf.gz"
        tmp2="${OUTDIR}/${tag}.tmp2.vcf.gz"

        bcftools view -v "$vtype" "$JOINT_VCF" -Oz -o "$tmp1"
        bcftools +setGT "$tmp1" -Oz -o "$tmp2" -- -t q -n . -i "FMT/DP<${dp}"
        bcftools view -m2 -M2 "$tmp2" -Oz -o "$out"
        bcftools index -t "$out"
        rm -f "$tmp1" "$tmp2"

        n=$(bcftools view -H "$out" | wc -l)
        echo "[INFO]   Sites in ${tag}: ${n}"
    done
done

echo "[STEP 4] Done."

# =============================================================================
# STEP 5: Count summary
# =============================================================================

echo ""
echo "[STEP 5] Generating count summary..."

module load python/3.10

python3 << 'PYEOF'
import os, csv, gzip
from collections import defaultdict

outdir = "/blue/baer/m.rifat/briggsae_ancestors/Joint_ancestor_call/summary"

tag_pairs = [
    ("snps.DP3",    "snp.DP3"),
    ("snps.DP5",    "snp.DP5"),
    ("snps.DP7",    "snp.DP7"),
    ("snps.DP10",   "snp.DP10"),
    ("indels.DP3",  "indel.DP3"),
    ("indels.DP5",  "indel.DP5"),
    ("indels.DP7",  "indel.DP7"),
    ("indels.DP10", "indel.DP10"),
]
out_tags = [t[1] for t in tag_pairs]

# ==========================================================================
# Parse sample name → canonical label
# ==========================================================================
def parse_sample(name):
    first = name.split('_')[0]
    anc = 'HK' if first.upper().startswith('HK') else \
          'PB' if first.upper().startswith('PB') else first

    if 'hifi' in name.lower() or 'G0' in name:
        return f"{anc}_LR"
    else:
        # SR DeepVariant gVCF — replicate from last token before .dv
        base = name.split('.dv')[0] if '.dv' in name else name
        rep  = 'R1' if base.split('_')[-1].endswith('1') else 'R2'
        return f"{anc}_SR_{rep}"

# ==========================================================================
# Load biallelic VCF → site_counts
# Strict filter: all 6 samples must be 0/0 or 1/1
# Indel filter:  |len(ALT) - len(REF)| <= 30
# ==========================================================================
def load_vcf(fpath, is_indel=False):
    site_counts = defaultdict(int)
    samples     = []
    n_total = n_skip_indel = n_skip_gt = 0

    opener = gzip.open if fpath.endswith('.gz') else open
    with opener(fpath, 'rt') as f:
        for line in f:
            if line.startswith('##'): continue
            if line.startswith('#CHROM'):
                raw = line.rstrip('\n').split('\t')[9:]
                samples = [parse_sample(s) for s in raw]
                continue

            n_total += 1
            fields = line.rstrip('\n').split('\t')
            ref, alt = fields[3], fields[4]

            # Indel size filter
            if is_indel and abs(len(alt) - len(ref)) > 30:
                n_skip_indel += 1
                continue

            fmt    = fields[8].split(':')
            gt_idx = fmt.index('GT') if 'GT' in fmt else 0

            gts = [fields[9+i].split(':')[gt_idx].replace('|','/')
                   for i in range(len(samples))]

            # Strict genotype filter
            def is_clean(gt):
                als = gt.split('/')
                return set(als) in ({'0'}, {'1'})

            if not all(is_clean(gt) for gt in gts):
                n_skip_gt += 1
                continue

            hom_set = frozenset(
                samples[i] for i, gt in enumerate(gts)
                if set(gt.split('/')) == {'1'}
            )
            if hom_set:
                site_counts[hom_set] += 1

    print(f"    Total sites      : {n_total}")
    if is_indel:
        print(f"    Skipped (>30bp)  : {n_skip_indel}")
    print(f"    Skipped (het/mis): {n_skip_gt}")
    print(f"    Sites kept       : {sum(site_counts.values())}")
    return site_counts, samples

# ==========================================================================
# Build row definitions (63 exact patterns + 6 grand totals)
# ==========================================================================
def build_row_defs():
    hk_lr = 'HK_LR'; pb_lr = 'PB_LR'
    hk_r1 = 'HK_SR_R1'; hk_r2 = 'HK_SR_R2'
    pb_r1 = 'PB_SR_R1'; pb_r2 = 'PB_SR_R2'
    S = [hk_lr, hk_r1, hk_r2, pb_lr, pb_r1, pb_r2]

    rows = []
    def sec(t):   rows.append((t, None))
    def row(l,p): rows.append((l, frozenset(p)))

    sec("=== SIZE 6: All 6 samples ===")
    row("All 6 (HK_LR+SR_R1+R2, PB_LR+SR_R1+R2)", S)

    sec("=== SIZE 5: 5 of 6 samples ===")
    sec("--- HK complete, PB missing one ---")
    row("HK all 3 + PB_LR + PB_SR_R1  (PB_SR_R2 absent)", [hk_lr,hk_r1,hk_r2,pb_lr,pb_r1])
    row("HK all 3 + PB_LR + PB_SR_R2  (PB_SR_R1 absent)", [hk_lr,hk_r1,hk_r2,pb_lr,pb_r2])
    row("HK all 3 + PB_SR_R1 + PB_SR_R2  (PB_LR absent)", [hk_lr,hk_r1,hk_r2,pb_r1,pb_r2])
    sec("--- PB complete, HK missing one ---")
    row("PB all 3 + HK_LR + HK_SR_R1  (HK_SR_R2 absent)", [pb_lr,pb_r1,pb_r2,hk_lr,hk_r1])
    row("PB all 3 + HK_LR + HK_SR_R2  (HK_SR_R1 absent)", [pb_lr,pb_r1,pb_r2,hk_lr,hk_r2])
    row("PB all 3 + HK_SR_R1 + HK_SR_R2  (HK_LR absent)", [pb_lr,pb_r1,pb_r2,hk_r1,hk_r2])

    sec("=== SIZE 4: 4 of 6 samples ===")
    sec("--- Both ancestors fully represented (2 from each) ---")
    row("HK_LR+SR_R1 + PB_LR+SR_R1  (HK_R2, PB_R2 absent)", [hk_lr,hk_r1,pb_lr,pb_r1])
    row("HK_LR+SR_R1 + PB_LR+SR_R2  (HK_R2, PB_R1 absent)", [hk_lr,hk_r1,pb_lr,pb_r2])
    row("HK_LR+SR_R1 + PB_SR_R1+R2  (HK_R2, PB_LR absent)", [hk_lr,hk_r1,pb_r1,pb_r2])
    row("HK_LR+SR_R2 + PB_LR+SR_R1  (HK_R1, PB_R2 absent)", [hk_lr,hk_r2,pb_lr,pb_r1])
    row("HK_LR+SR_R2 + PB_LR+SR_R2  (HK_R1, PB_R1 absent)", [hk_lr,hk_r2,pb_lr,pb_r2])
    row("HK_LR+SR_R2 + PB_SR_R1+R2  (HK_R1, PB_LR absent)", [hk_lr,hk_r2,pb_r1,pb_r2])
    row("HK_SR_R1+R2 + PB_LR+SR_R1  (HK_LR, PB_R2 absent)", [hk_r1,hk_r2,pb_lr,pb_r1])
    row("HK_SR_R1+R2 + PB_LR+SR_R2  (HK_LR, PB_R1 absent)", [hk_r1,hk_r2,pb_lr,pb_r2])
    row("HK_SR_R1+R2 + PB_SR_R1+R2  (both LR absent)",       [hk_r1,hk_r2,pb_r1,pb_r2])
    sec("--- All 3 from one ancestor + 1 from other ---")
    row("HK all 3 + PB_LR only",    [hk_lr,hk_r1,hk_r2,pb_lr])
    row("HK all 3 + PB_SR_R1 only", [hk_lr,hk_r1,hk_r2,pb_r1])
    row("HK all 3 + PB_SR_R2 only", [hk_lr,hk_r1,hk_r2,pb_r2])
    row("PB all 3 + HK_LR only",    [pb_lr,pb_r1,pb_r2,hk_lr])
    row("PB all 3 + HK_SR_R1 only", [pb_lr,pb_r1,pb_r2,hk_r1])
    row("PB all 3 + HK_SR_R2 only", [pb_lr,pb_r1,pb_r2,hk_r2])

    sec("=== SIZE 3: 3 of 6 samples ===")
    sec("--- All 3 from one ancestor only ---")
    row("HK all 3 only  (PB fully absent)", [hk_lr,hk_r1,hk_r2])
    row("PB all 3 only  (HK fully absent)", [pb_lr,pb_r1,pb_r2])
    sec("--- LR + one SR of same ancestor only ---")
    row("HK_LR + HK_SR_R1 only", [hk_lr,hk_r1])
    row("HK_LR + HK_SR_R2 only", [hk_lr,hk_r2])
    row("PB_LR + PB_SR_R1 only", [pb_lr,pb_r1])
    row("PB_LR + PB_SR_R2 only", [pb_lr,pb_r2])
    sec("--- Both SR of same ancestor only (no LR) ---")
    row("HK_SR_R1 + HK_SR_R2 only  (HK_LR absent)", [hk_r1,hk_r2])
    row("PB_SR_R1 + PB_SR_R2 only  (PB_LR absent)", [pb_r1,pb_r2])
    sec("--- Both LR + one SR ---")
    row("HK_LR + PB_LR + HK_SR_R1 only", [hk_lr,pb_lr,hk_r1])
    row("HK_LR + PB_LR + HK_SR_R2 only", [hk_lr,pb_lr,hk_r2])
    row("HK_LR + PB_LR + PB_SR_R1 only", [hk_lr,pb_lr,pb_r1])
    row("HK_LR + PB_LR + PB_SR_R2 only", [hk_lr,pb_lr,pb_r2])
    sec("--- One LR + two SR from different ancestors ---")
    row("HK_LR + HK_SR_R1 + PB_SR_R1 only", [hk_lr,hk_r1,pb_r1])
    row("HK_LR + HK_SR_R1 + PB_SR_R2 only", [hk_lr,hk_r1,pb_r2])
    row("HK_LR + HK_SR_R2 + PB_SR_R1 only", [hk_lr,hk_r2,pb_r1])
    row("HK_LR + HK_SR_R2 + PB_SR_R2 only", [hk_lr,hk_r2,pb_r2])
    row("PB_LR + PB_SR_R1 + HK_SR_R1 only", [pb_lr,pb_r1,hk_r1])
    row("PB_LR + PB_SR_R1 + HK_SR_R2 only", [pb_lr,pb_r1,hk_r2])
    row("PB_LR + PB_SR_R2 + HK_SR_R1 only", [pb_lr,pb_r2,hk_r1])
    row("PB_LR + PB_SR_R2 + HK_SR_R2 only", [pb_lr,pb_r2,hk_r2])
    sec("--- Cross-ancestor SR only (3 SR, no LR) ---")
    row("HK_SR_R1 + HK_SR_R2 + PB_SR_R1 only", [hk_r1,hk_r2,pb_r1])
    row("HK_SR_R1 + HK_SR_R2 + PB_SR_R2 only", [hk_r1,hk_r2,pb_r2])
    row("HK_SR_R1 + PB_SR_R1 + PB_SR_R2 only", [hk_r1,pb_r1,pb_r2])
    row("HK_SR_R2 + PB_SR_R1 + PB_SR_R2 only", [hk_r2,pb_r1,pb_r2])

    sec("=== SIZE 2: 2 of 6 samples ===")
    sec("--- Both LR only ---")
    row("HK_LR + PB_LR only  (no SR)", [hk_lr,pb_lr])
    sec("--- LR + SR cross ancestor ---")
    row("HK_LR + PB_SR_R1 only", [hk_lr,pb_r1])
    row("HK_LR + PB_SR_R2 only", [hk_lr,pb_r2])
    row("PB_LR + HK_SR_R1 only", [pb_lr,hk_r1])
    row("PB_LR + HK_SR_R2 only", [pb_lr,hk_r2])
    sec("--- SR cross ancestor (one each, no LR) ---")
    row("HK_SR_R1 + PB_SR_R1 only", [hk_r1,pb_r1])
    row("HK_SR_R1 + PB_SR_R2 only", [hk_r1,pb_r2])
    row("HK_SR_R2 + PB_SR_R1 only", [hk_r2,pb_r1])
    row("HK_SR_R2 + PB_SR_R2 only", [hk_r2,pb_r2])

    sec("=== SIZE 1: Single sample only ===")
    row("HK_LR only",    [hk_lr])
    row("HK_SR_R1 only", [hk_r1])
    row("HK_SR_R2 only", [hk_r2])
    row("PB_LR only",    [pb_lr])
    row("PB_SR_R1 only", [pb_r1])
    row("PB_SR_R2 only", [pb_r2])

    sec("=== GRAND TOTAL per sample (presence irrespective of others) ===")
    rows.append(("HK_LR total",    'grand', hk_lr))
    rows.append(("HK_SR_R1 total", 'grand', hk_r1))
    rows.append(("HK_SR_R2 total", 'grand', hk_r2))
    rows.append(("PB_LR total",    'grand', pb_lr))
    rows.append(("PB_SR_R1 total", 'grand', pb_r1))
    rows.append(("PB_SR_R2 total", 'grand', pb_r2))

    return rows

# ==========================================================================
# Main
# ==========================================================================
all_counts   = {}
row_defs     = build_row_defs()

for file_tag, out_tag in tag_pairs:
    vcf_path = os.path.join(outdir, f"{file_tag}.biallelic.vcf.gz")
    is_indel  = 'indel' in out_tag

    print(f"\n{'='*60}")
    print(f"Processing: {out_tag}")
    print(f"{'='*60}")

    if not os.path.exists(vcf_path):
        print(f"  ERROR: VCF not found: {vcf_path}")
        continue

    site_counts, _ = load_vcf(vcf_path, is_indel=is_indel)
    all_counts[out_tag] = site_counts

def exact(sc, pat): return sc.get(frozenset(pat), 0)
def grand(sc, lbl): return sum(v for k,v in sc.items() if lbl in k)

out_csv = os.path.join(outdir, "combined_count_summary_strict.csv")
with open(out_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['CATEGORY'] + out_tags)
    for r in row_defs:
        if r[1] is None:
            writer.writerow([r[0]] + ['']*len(out_tags))
        elif r[1] == 'grand':
            writer.writerow([r[0]] + [grand(all_counts.get(t,{}), r[2]) for t in out_tags])
        else:
            writer.writerow([r[0]] + [exact(all_counts.get(t,{}), r[1]) for t in out_tags])

print(f"\n[DONE] Written: {out_csv}", flush=True)
PYEOF

echo ""
echo "================================================================"
echo "[DONE] Full pipeline complete."
echo "Output: $OUTDIR/combined_count_summary_strict.csv"
echo "================================================================"
