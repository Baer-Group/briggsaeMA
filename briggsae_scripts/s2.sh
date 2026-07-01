#!/bin/sh
#SBATCH --job-name=s2 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=6 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=40gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=s2_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /orange/baer/briggsae/Baer_20251110/rtanny1_212714
### step 1: Decompress the fastq files:
gunzip *.fastq.gz
## Step 2: Quality trimming
module load fastp/
for i in `ls *_R1_001.fastq`; do
 file=$(basename $i "_R1_001.fastq")
 fastp --detect_adapter_for_pe \
 -i ${file}_R1_001.fastq \
 -I ${file}_R2_001.fastq \
 -o ${file}_R1_trimmed.fastq \
 -O ${file}_R2_trimmed.fastq
done
