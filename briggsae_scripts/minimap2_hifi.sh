#!/bin/sh
#SBATCH --job-name=minimap2 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=8 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=40gb # Memory limit
#SBATCH --time=480:00:00 # Time limit hrs:min:sec
#SBATCH --output=minimap2_%j.out # Standard output and error log
#SBATCH --account=juannanzhou --qos=juannanzhou


cd /orange/juannanzhou/Rifat_CB_raw/briggsae/Long_read_Fresh



module load minimap2

module load samtools



REF=20250626_c_briggsae_Feb2020.genome.fa



for fq in *filt.fastq.gz; do



    sample="${fq%.filt.fastq.gz}"



    minimap2 -t 8 -ax map-hifi "$REF" "$fq" | samtools sort -@ 8 -o "${sample}.bam"



    samtools index -@ 8 "${sample}.bam"



done
