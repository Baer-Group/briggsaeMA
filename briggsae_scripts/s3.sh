#!/bin/sh
#SBATCH --job-name=s3 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=6 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=40gb # Memory limit
#SBATCH --time=280:00:00 # Time limit hrs:min:sec
#SBATCH --output=s3_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /orange/baer/briggsae
### step 3: Bowtie2 alignment:
module load gcc/12.2.0/
module load bowtie2/2.4.5/
module load samtools/
module load bamtools/
#bowtie2-build 20250626_c_briggsae_Feb2020.genome.fa briggsae
#for i in `ls *_R1_trimmed.fastq`; do
#file=$(basename $i "_R1_trimmed.fastq")
#bowtie2 -x briggsae -1 ${file}_R1_trimmed.fastq -2 ${file}_R2_trimmed.fastq --phred33 -p 8 --very-sensitive-local -S ${file}.sam
#done
### step 4: Sam to Bam -- All alignments MQ<3/!PP are filtered out:
for i in `ls *.sam`
do
file=$(basename $i ".sam")
samtools view -b -q 3 -o ${file}.bam ${file}.sam
done
for i in *.bam
do
samtools sort -O bam -T tmp_ -o "${i/%.bam/.output.bam}" "$i"
done
for i in *.output.bam
do
bamtools filter -isProperPair true -in "$i" -out "${i/%.output.bam/.PP.bam}"
done
