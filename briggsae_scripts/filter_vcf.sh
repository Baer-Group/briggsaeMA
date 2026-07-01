#!/bin/bash
# =============================================================================
# VCF Filtering Script
# Filters: (1) heterozygous sites, (2) per-sample minority AF > 10%,
#          (3) sites with >=10 missing genotypes, (4) sites with no 1/1 calls
# =============================================================================

# --- Configuration ---
INPUT_DIR="."          # Directory containing input VCF files
OUTPUT_DIR="./filtered" # Output directory

mkdir -p "${OUTPUT_DIR}"

# --- Embedded Python filter script ---
PYTHON_SCRIPT=$(cat << 'EOF'
import sys
import os

def filter_vcf(input_file, output_file, min_missing=10, af_threshold=0.10):
    """
    Filters a VCF file based on the following criteria:
      1. Discard any site where at least one sample is heterozygous (0/1)
      2. Set a sample's genotype to missing (./.) if minority allele AF > 10%
           - For 0/0: minority AF = sum(ALT depths) / total depth
           - For 1/1: minority AF = REF depth / total depth
      3. Discard site if >= 10 samples have missing genotypes after step 2
      4. Discard site if no sample has a homozygous alt (1/1) genotype after step 2
    """
    sites_total     = 0
    sites_het       = 0
    sites_missing   = 0
    sites_no_var    = 0
    sites_kept      = 0

    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            # Pass header lines through unchanged
            if line.startswith('#'):
                fout.write(line)
                continue

            sites_total += 1
            fields = line.rstrip('\n').split('\t')

            if len(fields) < 9:
                fout.write(line)
                continue

            format_fields = fields[8].split(':')

            # Locate GT and AD indices in FORMAT
            try:
                gt_idx = format_fields.index('GT')
            except ValueError:
                sys.stderr.write(f"WARNING: GT field not found in FORMAT at line, skipping: {fields[0]}:{fields[1]}\n")
                continue
            try:
                ad_idx = format_fields.index('AD')
            except ValueError:
                sys.stderr.write(f"WARNING: AD field not found in FORMAT at line, skipping: {fields[0]}:{fields[1]}\n")
                continue

            samples = fields[9:]

            # ------------------------------------------------------------------
            # FILTER 1: Discard site if any sample is heterozygous
            # ------------------------------------------------------------------
            is_het = False
            for s in samples:
                sf = s.split(':')
                if gt_idx >= len(sf):
                    continue
                gt = sf[gt_idx].replace('|', '/')
                alleles = gt.split('/')
                # Heterozygous: not missing, not all same allele
                if '.' not in alleles and len(set(alleles)) > 1:
                    is_het = True
                    break
            if is_het:
                sites_het += 1
                continue

            # ------------------------------------------------------------------
            # FILTER 2: Per-sample minority AF check; set to ./. if AF > 10%
            # ------------------------------------------------------------------
            new_samples = []
            for s in samples:
                sf = s.split(':')
                if gt_idx >= len(sf):
                    new_samples.append(s)
                    continue

                gt = sf[gt_idx].replace('|', '/')
                alleles = gt.split('/')

                # Already missing — keep as is
                if '.' in alleles:
                    new_samples.append(s)
                    continue

                # Get AD values
                if ad_idx >= len(sf) or sf[ad_idx] in ('.', './.', ''):
                    new_samples.append(s)
                    continue

                try:
                    ad_vals = [int(x) for x in sf[ad_idx].split(',')]
                except ValueError:
                    new_samples.append(s)
                    continue

                total_depth = sum(ad_vals)
                if total_depth == 0:
                    sf[gt_idx] = './.'
                    new_samples.append(':'.join(sf))
                    continue

                allele_set = set(alleles)

                if allele_set == {'0'}:
                    # Homozygous REF: minority = all ALT depths combined
                    alt_depth = sum(ad_vals[1:]) if len(ad_vals) > 1 else 0
                    minority_af = alt_depth / total_depth

                elif allele_set == {'1'}:
                    # Homozygous ALT: minority = REF depth
                    ref_depth = ad_vals[0]
                    minority_af = ref_depth / total_depth

                else:
                    # Multi-allelic hom or unexpected — keep as is
                    new_samples.append(s)
                    continue

                if minority_af > af_threshold:
                    sf[gt_idx] = './.'
                    new_samples.append(':'.join(sf))
                else:
                    new_samples.append(s)

            # ------------------------------------------------------------------
            # FILTER 3: Discard site if >= 10 samples have missing genotypes
            # ------------------------------------------------------------------
            missing_count = 0
            for s in new_samples:
                sf = s.split(':')
                if gt_idx >= len(sf):
                    continue
                gt = sf[gt_idx].replace('|', '/')
                if '.' in gt.split('/'):
                    missing_count += 1

            if missing_count >= min_missing:
                sites_missing += 1
                continue

            # ------------------------------------------------------------------
            # FILTER 4: Discard site if no sample has a homozygous alt (1/1)
            # ------------------------------------------------------------------
            has_hom_alt = False
            for s in new_samples:
                sf = s.split(':')
                if gt_idx >= len(sf):
                    continue
                gt = sf[gt_idx].replace('|', '/').split('/')
                if '.' not in gt and len(set(gt)) == 1 and gt[0] != '0':
                    has_hom_alt = True
                    break

            if not has_hom_alt:
                sites_no_var += 1
                continue

            # ------------------------------------------------------------------
            # Site passed all filters — write out
            # ------------------------------------------------------------------
            fields[9:] = new_samples
            fout.write('\t'.join(fields) + '\n')
            sites_kept += 1

    # Print summary for this file
    print(f"  Total sites evaluated : {sites_total}")
    print(f"  Removed (heterozygous): {sites_het}")
    print(f"  Removed (>=10 missing): {sites_missing}")
    print(f"  Removed (no hom-alt)  : {sites_no_var}")
    print(f"  Sites kept            : {sites_kept}")

# Entry point when called from bash
if __name__ == '__main__':
    filter_vcf(sys.argv[1], sys.argv[2])
EOF
)

# --- Main loop over all VCF files ---
echo "========================================"
echo "Starting VCF filtering"
echo "Input  directory : ${INPUT_DIR}"
echo "Output directory : ${OUTPUT_DIR}"
echo "Filters applied  :"
echo "  1. Discard sites with any heterozygous genotype"
echo "  2. Set sample genotype to ./. if minority AF > 10%"
echo "  3. Discard sites with >= 10 missing genotypes"
echo "  4. Discard sites with no remaining 1/1 genotype"
echo "========================================"

for vcf in "${INPUT_DIR}"/*.biallelic.vcf; do
    [ -e "$vcf" ] || { echo "No VCF files found in ${INPUT_DIR}"; exit 1; }

    basename_vcf=$(basename "${vcf}" .vcf)
    output_vcf="${OUTPUT_DIR}/${basename_vcf}.af_filtered.vcf"

    echo ""
    echo "Processing: $(basename ${vcf})"
    echo "  -> Output: $(basename ${output_vcf})"

    python3 - "${vcf}" "${output_vcf}" << PYEOF
${PYTHON_SCRIPT}
PYEOF

    if [ $? -eq 0 ]; then
        echo "  Done."
    else
        echo "  ERROR processing ${vcf}" >&2
    fi
done

echo ""
echo "========================================"
echo "All files processed."
echo "Filtered VCFs written to: ${OUTPUT_DIR}"
echo "========================================"
