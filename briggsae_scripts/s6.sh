#!/bin/sh
#SBATCH --job-name=s6 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=9 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=60gb # Memory limit
#SBATCH --time=280:00:00 # Time limit hrs:min:sec
#SBATCH --output=s6_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /orange/baer/briggsae
module load gatk/
module load picard/
for i in `ls *.PP.RG.bam`; do
file=$(basename $i ".PP.RG.bam")
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar MarkDuplicates \
-I ${file}.PP.RG.bam \
-O ${file}.markduplicates.bam \
-M ${file}.markduplicates_metrics.txt \
--CREATE_INDEX true \
--REMOVE_DUPLICATES true \
--TMP_DIR test
done


