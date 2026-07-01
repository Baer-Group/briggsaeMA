#! /bin/sh
#SBATCH --job-name=gc # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=480:00:00 # Time limit hrs:min:sec
#SBATCH --output=gc_%j.out # Standard output and error log
#SBATCH --account=juannanzhou --qos=juannanzhou 
cd /orange/baer/briggsae
### step 5: Add read group information using picard:
module load samtools


echo -e "Sample\tMean_Coverage" > coverage_per_sample.tsv



for bam in *.markduplicates.bam

do

    sample=$(basename $bam .bam)



    cov=$(samtools coverage $bam | awk 'NR==2 {print $7}')



    echo -e "${sample}\t${cov}" >> coverage_per_sample.tsv

done
