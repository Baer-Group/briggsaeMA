#!/usr/bin/env bash
#SBATCH --job-name=dv # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=8 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=dv_%j.out # Standard output and error log
#SBATCH --account=juannanzhou --qos=juannanzhou-b
set -euo pipefail

cd /orange/juannanzhou/Rifat_CB_raw/briggsae/dummy_snp


# ── FIX 1: Use SLURM scratch instead of /tmp ─────────────────────────────────

# /tmp is shared and tiny on compute nodes; SLURM scratch is job-private & large

export TMPDIR="/scratch/local/${SLURM_JOB_ID}"

mkdir -p "$TMPDIR"

# ─────────────────────────────────────────────────────────────────────────────



module load deepvariant/1.10.0





REF="snp_dummy.fa"

THREADS=8



for BAM in *.bam; do

    base="${BAM%.bam}"

 # ── FIX 2: Skip if both output VCFs already exist ────────────────────────



    if [[ -f "${base}.dv.vcf.gz" && -f "${base}.dv.g.vcf.gz" ]]; then



        echo "[SKIP] ${base} — output files already exist, skipping."



        continue



    fi

    echo "Running DeepVariant on ${base}"



    deepvariant --model_type=PACBIO --ref="$REF" --reads="$BAM" --output_vcf="${base}.dv.vcf.gz" --output_gvcf="${base}.dv.g.vcf.gz" --num_shards="$THREADS"

done



