#!/bin/sh
#SBATCH --job-name=s9 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=6 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=50gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=s9_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /orange/baer/briggsae
### step 9: Call variants jointly using GenotypeGVCFs in GATK:
module load gatk/
gatk --java-options "-Xmx8g" GenotypeGVCFs -R 20250626_c_briggsae_Feb2020.genome.fa -V gendb://B2 -O B2.bp.GVCFs.vcf
