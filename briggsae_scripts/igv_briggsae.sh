#!/usr/bin/env bash
#SBATCH --job-name=igv_briggsae
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --time=96:00:00
#SBATCH --output=igv_briggsae_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou-b

set -euo pipefail

# =============================================================================
# IGV Automated Screenshot Pipeline — C. briggsae /orange
# =============================================================================
# Reads mutation_list CSV, maps each Sample to its B1/B2 BAM files,
# groups by HK (2xx) or PB (3xx), generates IGV batch scripts,
# and runs IGV headlessly via Xvfb to produce per-site screenshots.
#
# INPUT CSV columns: CHROM, POS, REF, ALT, Sample, B1, B2, TYPE
#
# OUTPUT STRUCTURE:
#   $SNAP_DIR/
#     HK/  {CHROM}_{POS}_{TYPE}_{SAMPLE}.png
#     PB/  {CHROM}_{POS}_{TYPE}_{SAMPLE}.png
#   $SHEET_DIR/
#     HK_contact.png
#     PB_contact.png
# =============================================================================

# =============================================================================
# SETTINGS — edit here only
# =============================================================================
WORKDIR="/orange/baer/briggsae"
INPUT_CSV="${WORKDIR}/mutation_list_5x_unique.csv"   # change to 10x etc. as needed
REF="${WORKDIR}/20250626_c_briggsae_Feb2020.genome.fa"
SNAP_DIR="${WORKDIR}/igv_screenshots"
BATCH_DIR="${WORKDIR}/igv_batches"
SHEET_DIR="${WORKDIR}/igv_sheets"
WINDOW=50       # bp to show on each side of the variant
IGV_MEM="20g"

mkdir -p "$SNAP_DIR" "$BATCH_DIR" "$SHEET_DIR"

# =============================================================================
# STEP 1: Parse CSV and generate per-group IGV batch files
# =============================================================================
echo "[INFO] Parsing CSV and generating IGV batch files..."

module load python/3.10
pip install openpyxl pandas Pillow --quiet --break-system-packages 2>/dev/null || true

python3 << PYEOF
import os, sys
import pandas as pd
from collections import defaultdict

INPUT_CSV = os.environ.get('INPUT_CSV', "${INPUT_CSV}")
WORKDIR   = os.environ.get('WORKDIR',   "${WORKDIR}")
SNAP_DIR  = os.environ.get('SNAP_DIR',  "${SNAP_DIR}")
BATCH_DIR = os.environ.get('BATCH_DIR', "${BATCH_DIR}")
REF       = os.environ.get('REF',       "${REF}")
WINDOW    = int(os.environ.get('WINDOW', "${WINDOW}"))

# ------------------------------------------------------------------
# Load CSV
# ------------------------------------------------------------------
df = pd.read_csv(INPUT_CSV)
df['Sample'] = df['Sample'].astype(str)
print(f"[INFO] Total sites in CSV: {len(df)}")

# ------------------------------------------------------------------
# Group assignment
# HK: sample starts with '2' or is an HK ancestor
# PB: sample starts with '3' or is a PB ancestor
# ------------------------------------------------------------------
HK_ANC_KEYS = {'HK104_ANC_1', 'HK104_ANC', 'HK_ANC'}
PB_ANC_KEYS = {'PB800_ANC__1', 'PB800_ANC', 'PB_ANC'}

def get_group(sample):
    if sample in HK_ANC_KEYS or sample.startswith('2'):
        return 'HK'
    elif sample in PB_ANC_KEYS or sample.startswith('3'):
        return 'PB'
    return None

# Collect all unique samples from CSV and classify
all_samples = df['Sample'].unique().tolist()
hk_samples  = [s for s in all_samples if get_group(s) == 'HK']
pb_samples  = [s for s in all_samples if get_group(s) == 'PB']
print(f"[INFO] HK samples ({len(hk_samples)}): {sorted(hk_samples)}")
print(f"[INFO] PB samples ({len(pb_samples)}): {sorted(pb_samples)}")

# ------------------------------------------------------------------
# BAM finder — B1 = L001, B2 = L002, no LR
# ------------------------------------------------------------------
def build_bam_cache(root):
    """Walk root once; map (sample_key, lane) -> bam_path."""
    cache = {}
    for dirpath, _, files in os.walk(root):
        for f in files:
            if not f.endswith('.bam') or 'markduplicates' not in f:
                continue
            full = os.path.join(dirpath, f)
            for lane in ('L001', 'L002'):
                if f'_{lane}' not in f:
                    continue
                idx = f.find('_POOLRET')
                if idx != -1:
                    key = f[:idx]
                    cache[(key, lane)] = full
    return cache

def find_bam(cache, sample, lane):
    # Direct lookup
    if (sample, lane) in cache:
        return cache[(sample, lane)]
    # ANC alias normalisation: strip trailing _1 / __1
    for alias in (sample.rstrip('_1').rstrip('_'), sample.replace('__1','').replace('_1','')):
        if (alias, lane) in cache:
            return cache[(alias, lane)]
    # Partial prefix match
    for (k_samp, k_lane), path in cache.items():
        if k_lane == lane and k_samp.startswith(sample):
            return path
    return None

print("[INFO] Building BAM cache (walking WORKDIR)...")
bam_cache = build_bam_cache(WORKDIR)
print(f"[INFO] BAMs indexed: {len(bam_cache)}")

# Report BAM resolution for all samples
def get_sample_bams(sample):
    bams = []
    for lane in ('L001', 'L002'):
        b = find_bam(bam_cache, sample, lane)
        status = os.path.basename(b) if b else 'NOT FOUND'
        print(f"  {sample}:{lane} -> {status}")
        if b:
            bams.append(b)
    return bams

print("\n[INFO] HK BAM resolution:")
hk_bams_all = []
for s in sorted(hk_samples):
    hk_bams_all.extend(get_sample_bams(s))
hk_bams_all = sorted(set(hk_bams_all))

print("\n[INFO] PB BAM resolution:")
pb_bams_all = []
for s in sorted(pb_samples):
    pb_bams_all.extend(get_sample_bams(s))
pb_bams_all = sorted(set(pb_bams_all))

print(f"\n[INFO] HK group: {len(hk_bams_all)} BAMs total")
print(f"[INFO] PB group: {len(pb_bams_all)} BAMs total")

# ------------------------------------------------------------------
# Group sites from CSV
# ------------------------------------------------------------------
groups = defaultdict(list)

for _, row in df.iterrows():
    sample = str(row['Sample'])
    chrom  = str(row['CHROM'])
    pos    = int(row['POS'])
    vtype  = str(row['TYPE'])
    b1_val = str(row['B1'])
    b2_val = str(row['B2'])

    grp = get_group(sample)
    if grp is None:
        print(f"  [SKIP] Unknown group for sample {sample} at {chrom}:{pos}")
        continue

    bams     = hk_bams_all if grp == 'HK' else pb_bams_all
    snap_sub = os.path.join(SNAP_DIR, grp)
    os.makedirs(snap_sub, exist_ok=True)

    # Label encodes sample + replicate status for easy reading
    b1_tag = 'B1mut' if b1_val not in ('absent','0/0','missing','') else 'B1ref'
    b2_tag = 'B2mut' if b2_val not in ('absent','0/0','missing','') else 'B2ref'
    label  = f"{chrom}_{pos}_{vtype}_{sample}_{b1_tag}_{b2_tag}"

    groups[grp].append({
        'chrom':    chrom,
        'pos':      pos,
        'label':    label,
        'snap_sub': snap_sub,
        'bams':     bams,
    })

print(f"\n[INFO] HK sites to screenshot: {len(groups['HK'])}")
print(f"[INFO] PB sites to screenshot: {len(groups['PB'])}")

# ------------------------------------------------------------------
# Write one IGV batch file per group
# ------------------------------------------------------------------
batch_files = []
for grp, sites in groups.items():
    if not sites:
        continue
    batch_path = os.path.join(BATCH_DIR, f"batch_{grp}.igv")
    bam_list   = sites[0]['bams']

    with open(batch_path, 'w') as f:
        f.write("new\n")
        f.write(f"genome {REF}\n")
        for bam in bam_list:
            f.write(f"load {bam}\n")
        f.write("preference SAM.SHOW_SOFT_CLIPPED true\n")
        f.write("preference SAM.COLOR_BY READ_STRAND\n")
        f.write("preference SAM.SHOW_ALL_BASES true\n")
        for site in sorted(sites, key=lambda x: (x['chrom'], x['pos'])):
            start = max(1, site['pos'] - WINDOW)
            end   = site['pos'] + WINDOW
            f.write(f"snapshotDirectory {site['snap_sub']}\n")
            f.write(f"goto {site['chrom']}:{start}-{end}\n")
            f.write("sort position\n")
            f.write("squish\n")
            f.write(f"snapshot {site['label']}.png\n")
        f.write("exit\n")

    batch_files.append(batch_path)
    print(f"[INFO] Batch written: {batch_path} ({len(sites)} sites)")

# Write batch list
list_path = os.path.join(BATCH_DIR, "batch_list.txt")
with open(list_path, 'w') as f:
    for b in batch_files:
        f.write(b + '\n')
print(f"[INFO] {len(batch_files)} batch files ready (HK + PB)")
PYEOF

# =============================================================================
# STEP 2: Run IGV headlessly via Xvfb for each batch file
# =============================================================================
echo "[INFO] Starting IGV screenshot runs..."

module load igv/2.18.4
IGV_CMD="igv.sh"
echo "[INFO] IGV 2.18.4 loaded"

# Start virtual display
DISPLAY_NUM=99
Xvfb :${DISPLAY_NUM} -screen 0 1920x1080x24 &
XVFB_PID=$!
export DISPLAY=:${DISPLAY_NUM}
sleep 3
echo "[INFO] Virtual display started (PID: $XVFB_PID)"

BATCH_LIST="${BATCH_DIR}/batch_list.txt"
TOTAL=$(wc -l < "$BATCH_LIST")
COUNT=0

while IFS= read -r batch_file; do
    COUNT=$((COUNT + 1))
    echo "[INFO] Running batch $COUNT/$TOTAL: $(basename $batch_file)"
    ${IGV_CMD} \
        --batch "$batch_file" \
    && echo "  [OK] Batch $COUNT complete" \
    || echo "  [WARN] Batch $COUNT exited with error — continuing"
done < "$BATCH_LIST"

kill $XVFB_PID 2>/dev/null || true
echo "[INFO] Virtual display stopped"

# =============================================================================
# STEP 3: Generate contact sheets (grid of screenshots per group)
# =============================================================================
echo "[INFO] Generating contact sheets..."
module load python/3.10

python3 << PYEOF3
import os
from PIL import Image, ImageDraw
import math

snap_dir  = "${SNAP_DIR}"
sheet_dir = "${SHEET_DIR}"
cols      = 5   # screenshots per row in contact sheet

for subdir in sorted(os.listdir(snap_dir)):
    subpath = os.path.join(snap_dir, subdir)
    if not os.path.isdir(subpath):
        continue

    imgs = sorted([f for f in os.listdir(subpath) if f.endswith('.png')])
    if not imgs:
        print(f"  [SKIP] {subdir} — no screenshots")
        continue

    sample_img = Image.open(os.path.join(subpath, imgs[0]))
    w, h = sample_img.size
    label_h = 30

    rows    = math.ceil(len(imgs) / cols)
    sheet_w = cols * w
    sheet_h = rows * (h + label_h)

    sheet = Image.new('RGB', (sheet_w, sheet_h), color='white')
    draw  = ImageDraw.Draw(sheet)

    for idx, img_file in enumerate(imgs):
        img = Image.open(os.path.join(subpath, img_file))
        r, c = divmod(idx, cols)
        x = c * w
        y = r * (h + label_h)
        sheet.paste(img, (x, y))
        # Label: CHROM_POS_TYPE from filename
        parts = img_file.replace('.png', '').split('_')
        label = '_'.join(parts[:4]) if len(parts) >= 4 else img_file[:35]
        draw.text((x + 4, y + h + 2), label, fill='black')

    out_path = os.path.join(sheet_dir, f"{subdir}_contact.png")
    sheet.save(out_path)
    print(f"  [OK] {subdir}: {len(imgs)} sites -> {os.path.basename(out_path)}")

print("[INFO] All contact sheets generated")
PYEOF3

# =============================================================================
# STEP 4: Summary
# =============================================================================
N_SNAPS=$(find "${SNAP_DIR}" -name "*.png" | wc -l)
echo ""
echo "================================================================"
echo " IGV Screenshot Pipeline Complete"
echo "================================================================"
echo " Input CSV   : ${INPUT_CSV}"
echo " Screenshots : ${SNAP_DIR}/"
echo " Contact sheets: ${SHEET_DIR}/"
echo " Batch files : ${BATCH_DIR}/"
echo " Total screenshots produced: ${N_SNAPS}"
echo "================================================================"
