#!/usr/bin/env bash
# =============================================================================
#  af_fpfn_v6.sh — BAM-level allele-frequency audit for C. briggsae MA lines
#
#  Goal: For every mutation in the pipeline's call list, go back to the raw
#  BAM reads and check whether the allele frequency at that site is consistent
#  with the call.  Each site has two replicates (B1 = L001, B2 = L002).
#
#  Verdict per replicate:
#    TRUE_POS  pipeline called mutation  AND  BAM AF ≥ 0.80, DP ≥ 5
#    FP        pipeline called mutation  BUT  BAM AF < 0.80 or DP < 5
#    FN        pipeline did NOT call     AND  BAM AF ≥ 0.80, DP ≥ 5
#    TRUE_NEG  pipeline did NOT call     AND  BAM AF < 0.80 (or low depth)
#    NO_BAM    BAM or index file missing
#    NO_DATA   pileup found no usable reads (dp == 0)
#
#  BEFORE FIRST RUN — pre-install packages on the login node once:
#    module load python/3.10
#    pip install pysam openpyxl pandas --user --quiet
#  Compute nodes inherit ~/.local/lib/python3.10/site-packages automatically.
# =============================================================================
#SBATCH --job-name=af_fpfn_v6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --output=af_fpfn_v6_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou

set -euo pipefail

module load python/3.10

WORKDIR="/orange/baer/briggsae"
INPUT_CSV="${WORKDIR}/mutation_list_5x_unique.csv"
OUTDIR="${WORKDIR}/af_fpfn_output"
mkdir -p "$OUTDIR"

# Verify packages are available before committing to the full job
python3 -c "import pysam, openpyxl, pandas" 2>/dev/null || {
    echo "[ERROR] Required Python packages not found."
    echo "        Run once on the login node:"
    echo "          module load python/3.10"
    echo "          pip install pysam openpyxl pandas --user --quiet"
    exit 1
}

python3 << 'PYEOF'
import os
import sys
from collections import Counter, defaultdict

import pysam
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

# =============================================================================
# CONFIG
# =============================================================================
WORKDIR   = "/orange/baer/briggsae"
INPUT_CSV = os.path.join(WORKDIR, "mutation_list_5x_unique.csv")
OUTDIR    = os.path.join(WORKDIR, "af_fpfn_output")
os.makedirs(OUTDIR, exist_ok=True)

MIN_DP = 5      # minimum depth to call a site informative
MIN_AF = 0.80   # minimum ALT allele frequency to call a site mutant


# =============================================================================
# BAM CACHE
#   Walk WORKDIR once.  Build {(sample_key, "L001"|"L002"): bam_path}.
#   Filename pattern: {SAMPLE}_POOLRET91_S{XX}_L{001|002}.markduplicates.bam
# =============================================================================
def build_bam_cache(root):
    cache      = {}
    duplicates = defaultdict(list)

    for dirpath, _, files in os.walk(root):
        for fname in files:
            if not fname.endswith(".bam"):
                continue
            if "markduplicates" not in fname:
                continue

            full = os.path.join(dirpath, fname)

            lane = None
            if "_L001" in fname:
                lane = "L001"
            elif "_L002" in fname:
                lane = "L002"
            else:
                continue

            idx = fname.find("_POOLRET")
            if idx == -1:
                continue

            sample_key = fname[:idx].strip()
            key = (sample_key, lane)

            if key in cache:
                duplicates[key].append(full)
            else:
                cache[key] = full

    if duplicates:
        print("[WARN] Duplicate BAM keys (keeping first found):")
        for key, extras in duplicates.items():
            print(f"  {key}: primary={cache[key]}")
            for e in extras:
                print(f"         extra  ={e}")

    return cache


def _normalize_sample(s):
    """Strip trailing _1 / __1 that appear on ANC sample names in the CSV."""
    if s.endswith("__1"):
        return s[:-3]
    if s.endswith("_1"):
        return s[:-2]
    return s


def find_bam(cache, sample, lane):
    """
    Exact match → normalized-suffix match.
    No prefix/partial matching (would silently resolve the wrong sample).
    """
    s = str(sample).strip()
    if (s, lane) in cache:
        return cache[(s, lane)]
    s2 = _normalize_sample(s)
    if (s2, lane) in cache:
        return cache[(s2, lane)]
    return None


# =============================================================================
# VARIANT TYPE
# =============================================================================
def variant_type(ref_al, alt_al):
    if len(ref_al) == 1 and len(alt_al) == 1:
        return "SNP"
    if len(alt_al) > len(ref_al):
        return "INS"
    if len(ref_al) > len(alt_al):
        return "DEL"
    return "COMPLEX"


# =============================================================================
# SOFT-CLIP HELPER  (DEL branch only)
#   Walk the CIGAR of an AlignmentSegment and return True when a soft-clip
#   (op=4) begins immediately after query position qpos.
#   This is the signature of a long deletion: the aligner anchors at qpos
#   and soft-clips the downstream bases rather than writing a D operation.
# =============================================================================
CONSUMES_QUERY = {0, 1, 4, 7, 8}   # M I S = X

def _softclip_after_anchor(aln, qpos):
    q = 0
    cigar = aln.cigartuples
    if not cigar:
        return False
    for op, length in cigar:
        if op == 4 and q == qpos + 1:  # soft-clip starts right after anchor
            return True
        if op in CONSUMES_QUERY:
            q += length
        if q > qpos + 1:
            break
    return False


# =============================================================================
# CORE COUNTING FUNCTION
#
#  Returns a dict:
#    ref_ct  : int | "."    reads supporting REF
#    alt_ct  : int | "."    reads supporting ALT
#    dp      : int | "."    ref_ct + alt_ct  (only informative reads)
#    af      : float | None alt_ct / dp;  None when dp == 0 or BAM unavailable
#    note    : str          "OK" or diagnostic flag(s)
#
#  Allele detection per type
#  ─────────────────────────
#  SNP:
#    pileup at anchor; base == ALT → alt, base == REF → ref.
#
#  INS (e.g. REF=A ALT=AAT → ins_seq="AT"):
#    fetch() all reads overlapping anchor.
#    For each, find query position (qpos) at anchor via get_aligned_pairs().
#    Check read.query_sequence[qpos+1 : qpos+1+ins_len] == ins_seq.
#    This single check works for both:
#      • CIGAR I reads  — inserted bases are at qpos+1 in the query sequence
#      • Soft-clipped reads — same position; aligner just called them S not I
#    No window guessing; exact positional check only.
#
#  DEL (e.g. REF=GTT ALT=G → del_len=2):
#    pileup at anchor.
#    ① read.indel == -del_len  → exact CIGAR D  → alt
#    ② soft-clip after anchor  → long DEL  → alt
#    ③ read.indel == 0, no soft-clip → ref
#    Reads ending at anchor are skipped (ambiguous).
# =============================================================================
def get_counts(bam_path, chrom, pos, ref_al, alt_al):
    empty = {"ref_ct": ".", "alt_ct": ".", "dp": ".", "af": None, "note": ""}

    # ── BAM / index checks ────────────────────────────────────────────────
    if not bam_path:
        return {**empty, "note": "NO_BAM"}
    if not os.path.exists(bam_path):
        return {**empty, "note": "NO_BAM"}
    if not (os.path.exists(bam_path + ".bai") or
            os.path.exists(bam_path.replace(".bam", ".bai"))):
        return {**empty, "note": "NO_INDEX"}

    vtype    = variant_type(ref_al, alt_al)
    anchor_0 = pos - 1          # 0-based
    ref_base = ref_al[0].upper()
    alt_base = alt_al[0].upper()

    ref_ct = 0
    alt_ct = 0
    notes  = []

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as exc:
        return {**empty, "note": f"BAM_OPEN_ERROR:{exc}"}

    try:
        with bam:

            # ── INS ──────────────────────────────────────────────────────
            if vtype == "INS":
                ins_seq = alt_al[len(ref_al):].upper()   # e.g. "AT" for REF=A ALT=AAT
                ins_len = len(ins_seq)
                n_reads = 0

                for read in bam.fetch(chrom, anchor_0, anchor_0 + 1):
                    if (read.is_unmapped or read.is_duplicate
                            or read.is_secondary or read.is_supplementary
                            or read.is_qcfail):
                        continue
                    if read.query_sequence is None or read.cigartuples is None:
                        continue
                    if (read.reference_start > anchor_0 or
                            read.reference_end is None or
                            read.reference_end <= anchor_0):
                        continue

                    # Find query index of the anchor base
                    qpos = None
                    for qp, rp in read.get_aligned_pairs(matches_only=False):
                        if rp == anchor_0 and qp is not None:
                            qpos = qp
                            break

                    if qpos is None:
                        continue
                    # Need at least ins_len bases after anchor in query
                    if qpos + 1 + ins_len > len(read.query_sequence):
                        continue
                    # Anchor must be REF base
                    if read.query_sequence[qpos].upper() != ref_base:
                        continue

                    n_reads += 1

                    # Exact positional check — works for CIGAR I AND soft-clip
                    candidate = read.query_sequence[qpos + 1: qpos + 1 + ins_len].upper()
                    if candidate == ins_seq:
                        alt_ct += 1
                    else:
                        ref_ct += 1

                if n_reads == 0:
                    notes.append("NO_READS")

            # ── DEL ──────────────────────────────────────────────────────
            elif vtype == "DEL":
                del_len   = len(ref_al) - len(alt_al)
                found_col = False

                for pcol in bam.pileup(chrom, anchor_0, anchor_0 + 1,
                                       truncate=True,
                                       stepper="samtools",
                                       min_base_quality=0,
                                       min_mapping_quality=0,
                                       ignore_overlaps=False,
                                       ignore_orphans=False):
                    if pcol.reference_pos != anchor_0:
                        continue
                    found_col = True

                    for pread in pcol.pileups:
                        aln = pread.alignment

                        if aln.is_secondary or aln.is_supplementary or aln.is_qcfail:
                            continue
                        if pread.is_refskip:
                            continue
                        if aln.query_sequence is None:
                            continue

                        # ① Exact-length CIGAR D
                        if pread.indel == -del_len:
                            alt_ct += 1
                            continue

                        # Reads with is_del at anchor have a deletion of a
                        # *different* length (multi-allelic or complex) → skip
                        if pread.is_del:
                            continue

                        qpos = pread.query_position
                        if qpos is None:
                            continue
                        if aln.query_sequence[qpos].upper() != ref_base:
                            continue

                        # Read ends at anchor — can't distinguish REF from DEL
                        if qpos + 1 >= aln.query_length:
                            continue

                        # ② Soft-clip immediately after anchor → long DEL
                        if aln.cigartuples and _softclip_after_anchor(aln, qpos):
                            alt_ct += 1
                            continue

                        # ③ No deletion evidence → REF
                        if pread.indel == 0:
                            ref_ct += 1

                if not found_col:
                    notes.append("NO_PILEUP_COL")

            # ── SNP ──────────────────────────────────────────────────────
            elif vtype == "SNP":
                found_col = False

                for pcol in bam.pileup(chrom, anchor_0, anchor_0 + 1,
                                       truncate=True,
                                       stepper="samtools",
                                       min_base_quality=0,
                                       min_mapping_quality=0,
                                       ignore_overlaps=False,
                                       ignore_orphans=False):
                    if pcol.reference_pos != anchor_0:
                        continue
                    found_col = True

                    for pread in pcol.pileups:
                        aln = pread.alignment

                        if aln.is_secondary or aln.is_supplementary or aln.is_qcfail:
                            continue
                        if pread.is_refskip or pread.is_del:
                            continue

                        qpos = pread.query_position
                        if qpos is None or aln.query_sequence is None:
                            continue

                        base = aln.query_sequence[qpos].upper()
                        if base == alt_base:
                            alt_ct += 1
                        elif base == ref_base:
                            ref_ct += 1
                        # other bases: not counted (sequencing error at anchor)

                if not found_col:
                    notes.append("NO_PILEUP_COL")

            # ── COMPLEX ──────────────────────────────────────────────────
            else:
                notes.append("COMPLEX_SKIPPED")

    except Exception as exc:
        return {**empty, "note": f"PILEUP_ERROR:{exc}"}

    dp = ref_ct + alt_ct

    # af is None (not 0.0) when dp == 0 so verdict() correctly returns NO_DATA
    af = (alt_ct / dp) if dp > 0 else None

    if dp < MIN_DP:
        notes.append("LOW_DP")

    note_str = ";".join(notes) if notes else "OK"
    return {"ref_ct": ref_ct, "alt_ct": alt_ct, "dp": dp, "af": af, "note": note_str}


# =============================================================================
# VERDICT
# =============================================================================
_NOT_CALLED = {"", ".", "na", "nan", "none", "absent", "missing",
               "0/0", "0|0", "ref", "reference"}

def is_called(field):
    return str(field).strip().lower() not in _NOT_CALLED

def verdict(b_col, counts):
    note = counts["note"]
    af   = counts["af"]
    dp   = counts["dp"]

    if note.startswith("NO_BAM"):
        return "NO_BAM"
    if note.startswith("NO_INDEX"):
        return "NO_INDEX"
    if note.startswith("BAM_OPEN_ERROR"):
        return "NO_BAM"
    if note.startswith("PILEUP_ERROR"):
        return "PILEUP_ERROR"
    if af is None:                  # dp == 0 or COMPLEX_SKIPPED
        return "NO_DATA"

    called = is_called(b_col)
    passes = (dp >= MIN_DP) and (af >= MIN_AF)

    if called and passes:      return "TRUE_POS"
    if called and not passes:  return "FP"
    if not called and passes:  return "FN"
    return "TRUE_NEG"


def fmt(x):
    """Format af to 3 dp, or '.' if None."""
    return f"{x:.3f}" if x is not None else "."


# =============================================================================
# LOAD INPUT
# =============================================================================
print("[INFO] Reading input CSV...")
df = pd.read_csv(INPUT_CSV)
print(f"[INFO] {len(df)} rows")

for col in ["CHROM", "POS", "REF", "ALT", "Sample", "B1", "B2"]:
    if col not in df.columns:
        sys.exit(f"[ERROR] Required column missing from CSV: {col}")

if "TYPE" not in df.columns:
    df["TYPE"] = [variant_type(str(r.REF), str(r.ALT)) for _, r in df.iterrows()]

# Coerce to string so NaN B1/B2 fields don't crash is_called()
df["B1"] = df["B1"].fillna("").astype(str)
df["B2"] = df["B2"].fillna("").astype(str)

print("[INFO] Building BAM cache...")
bam_cache = build_bam_cache(WORKDIR)
print(f"[INFO] {len(bam_cache)} BAMs indexed")

print("[INFO] Sample → BAM resolution:")
for s in sorted(df["Sample"].astype(str).unique()):
    for lane in ("L001", "L002"):
        p = find_bam(bam_cache, s, lane)
        print(f"  {s} [{lane}] → {os.path.basename(p) if p else 'NOT FOUND'}")


# =============================================================================
# EXCEL SETUP
# =============================================================================
HEADER_FILL = PatternFill("solid", fgColor="2E4057")
HEADER_FONT = Font(bold=True, color="FFFFFF", name="Arial")
BODY_FONT   = Font(name="Arial")

FILL = {
    "TRUE_POS"    : PatternFill("solid", fgColor="C8E6C9"),   # green
    "FP"          : PatternFill("solid", fgColor="FFB3B3"),   # red
    "FN"          : PatternFill("solid", fgColor="FFE0B2"),   # orange
    "TRUE_NEG"    : PatternFill("solid", fgColor="F5F5F5"),   # light grey
    "NO_BAM"      : PatternFill("solid", fgColor="D9D9D9"),   # grey
    "NO_INDEX"    : PatternFill("solid", fgColor="D9D9D9"),
    "NO_DATA"     : PatternFill("solid", fgColor="FFF2CC"),   # yellow
    "PILEUP_ERROR": PatternFill("solid", fgColor="FFF2CC"),
}
NOTE_WARN = PatternFill("solid", fgColor="FFF2CC")

wb = Workbook()
ws = wb.active
ws.title = "AF_FP_FN"

COLS = [
    "CHROM", "POS", "REF", "ALT", "Sample", "TYPE",
    "B1_called", "B1_REF", "B1_ALT", "B1_DP", "B1_AF", "B1_note", "B1_verdict",
    "B2_called", "B2_REF", "B2_ALT", "B2_DP", "B2_AF", "B2_note", "B2_verdict",
]

# Derive verdict/note column indices once — no magic numbers later
B1_VERD_COL = COLS.index("B1_verdict") + 1   # 1-based
B2_VERD_COL = COLS.index("B2_verdict") + 1
B1_NOTE_COL = COLS.index("B1_note")    + 1
B2_NOTE_COL = COLS.index("B2_note")    + 1

for ci, name in enumerate(COLS, 1):
    cell            = ws.cell(row=1, column=ci, value=name)
    cell.font       = HEADER_FONT
    cell.fill       = HEADER_FILL
    cell.alignment  = Alignment(horizontal="center")

COL_WIDTHS = {
    "CHROM": 10, "POS": 12, "REF": 20, "ALT": 30, "Sample": 14, "TYPE": 8,
    "B1_called": 16, "B1_REF": 7, "B1_ALT": 7, "B1_DP": 7, "B1_AF": 8,
    "B1_note": 22, "B1_verdict": 12,
    "B2_called": 16, "B2_REF": 7, "B2_ALT": 7, "B2_DP": 7, "B2_AF": 8,
    "B2_note": 22, "B2_verdict": 12,
}
for ci, name in enumerate(COLS, 1):
    ws.column_dimensions[get_column_letter(ci)].width = COL_WIDTHS.get(name, 10)

ws.freeze_panes = "G2"


# =============================================================================
# MAIN LOOP
# =============================================================================
total    = len(df)
sum_b1   = Counter()
sum_b2   = Counter()
type_b1  = defaultdict(Counter)
type_b2  = defaultdict(Counter)

for ri, (_, row) in enumerate(df.iterrows(), start=2):
    chrom  = str(row["CHROM"]).strip()
    pos    = int(row["POS"])
    ref_al = str(row["REF"]).strip().upper()
    alt_al = str(row["ALT"]).strip().upper()
    sample = str(row["Sample"]).strip()
    vtype  = str(row["TYPE"]).strip()
    b1_col = str(row["B1"]).strip()
    b2_col = str(row["B2"]).strip()

    if (ri - 2) % 50 == 0:
        print(f"[INFO] {ri - 1}/{total}  {chrom}:{pos}  {sample}  {vtype}")

    bam1 = find_bam(bam_cache, sample, "L001")
    bam2 = find_bam(bam_cache, sample, "L002")

    c1 = get_counts(bam1, chrom, pos, ref_al, alt_al)
    c2 = get_counts(bam2, chrom, pos, ref_al, alt_al)

    v1 = verdict(b1_col, c1)
    v2 = verdict(b2_col, c2)

    sum_b1[v1] += 1;  sum_b2[v2] += 1
    type_b1[vtype][v1] += 1;  type_b2[vtype][v2] += 1

    row_vals = [
        chrom, pos, ref_al, alt_al, sample, vtype,
        b1_col, c1["ref_ct"], c1["alt_ct"], c1["dp"], fmt(c1["af"]), c1["note"], v1,
        b2_col, c2["ref_ct"], c2["alt_ct"], c2["dp"], fmt(c2["af"]), c2["note"], v2,
    ]

    for ci, val in enumerate(row_vals, 1):
        cell      = ws.cell(row=ri, column=ci, value=val)
        cell.font = BODY_FONT
        if ci == B1_VERD_COL:
            cell.fill = FILL.get(v1, FILL["TRUE_NEG"])
        elif ci == B2_VERD_COL:
            cell.fill = FILL.get(v2, FILL["TRUE_NEG"])
        elif ci in (B1_NOTE_COL, B2_NOTE_COL) and str(val) not in ("OK", ""):
            cell.fill = NOTE_WARN


# =============================================================================
# SUMMARY SHEET
# =============================================================================
ws2 = wb.create_sheet("SUMMARY")

def hdr(ws, r, c, v):
    cell = ws.cell(row=r, column=c, value=v)
    cell.font = HEADER_FONT
    cell.fill = HEADER_FILL
    cell.alignment = Alignment(horizontal="center")

row = 1
# ── Overall counts ────────────────────────────────────────────────────────
hdr(ws2, row, 1, "Lane"); hdr(ws2, row, 2, "Verdict"); hdr(ws2, row, 3, "Count")
row += 1
for v, n in sorted(sum_b1.items()):
    ws2.cell(row=row, column=1, value="B1")
    ws2.cell(row=row, column=2, value=v)
    ws2.cell(row=row, column=3, value=n)
    row += 1
for v, n in sorted(sum_b2.items()):
    ws2.cell(row=row, column=1, value="B2")
    ws2.cell(row=row, column=2, value=v)
    ws2.cell(row=row, column=3, value=n)
    row += 1

# ── Per-TYPE breakdown ────────────────────────────────────────────────────
row += 1
hdr(ws2, row, 1, "Lane"); hdr(ws2, row, 2, "TYPE")
hdr(ws2, row, 3, "Verdict"); hdr(ws2, row, 4, "Count")
row += 1
all_types    = sorted(set(type_b1) | set(type_b2))
all_verdicts = sorted({"TRUE_POS", "FP", "FN", "TRUE_NEG", "NO_DATA",
                        "NO_BAM", "NO_INDEX", "PILEUP_ERROR"})
for vtype_key in all_types:
    for vname in all_verdicts:
        b1n = type_b1[vtype_key].get(vname, 0)
        b2n = type_b2[vtype_key].get(vname, 0)
        if b1n > 0:
            ws2.cell(row=row, column=1, value="B1")
            ws2.cell(row=row, column=2, value=vtype_key)
            ws2.cell(row=row, column=3, value=vname)
            ws2.cell(row=row, column=4, value=b1n)
            row += 1
        if b2n > 0:
            ws2.cell(row=row, column=1, value="B2")
            ws2.cell(row=row, column=2, value=vtype_key)
            ws2.cell(row=row, column=3, value=vname)
            ws2.cell(row=row, column=4, value=b2n)
            row += 1

# ── Run parameters ────────────────────────────────────────────────────────
row += 1
hdr(ws2, row, 1, "Parameter"); hdr(ws2, row, 2, "Value")
row += 1
for k, v in [("MIN_DP", MIN_DP), ("MIN_AF", MIN_AF),
              ("INPUT_CSV", INPUT_CSV), ("Script", "af_fpfn_v6")]:
    ws2.cell(row=row, column=1, value=k)
    ws2.cell(row=row, column=2, value=str(v))
    row += 1

for col_letter, w in zip("ABCD", [16, 20, 14, 10]):
    ws2.column_dimensions[col_letter].width = w


# =============================================================================
# SAVE + PRINT SUMMARY
# =============================================================================
out_path = os.path.join(OUTDIR, "af_fpfn_results_v6.xlsx")
wb.save(out_path)
print(f"\n[DONE] {out_path}")

print("\n=== Verdict summary ===")
print(f"{'Verdict':<14}  {'B1':>6}  {'B2':>6}")
print("-" * 32)
for v in all_verdicts:
    b1n = sum_b1.get(v, 0)
    b2n = sum_b2.get(v, 0)
    if b1n or b2n:
        print(f"{v:<14}  {b1n:>6}  {b2n:>6}")

print("\n=== Per-TYPE verdict summary ===")
for vtype_key in all_types:
    print(f"\n  {vtype_key}:")
    for vname in all_verdicts:
        b1n = type_b1[vtype_key].get(vname, 0)
        b2n = type_b2[vtype_key].get(vname, 0)
        if b1n or b2n:
            print(f"    {vname:<14}  B1={b1n:>4}  B2={b2n:>4}")

PYEOF

echo "================================================================"
echo " Done: ${OUTDIR}/af_fpfn_results_v6.xlsx"
echo "================================================================"
