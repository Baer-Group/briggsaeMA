#! /bin/sh
#SBATCH --job-name=cov # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=2 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=480:00:00 # Time limit hrs:min:sec
#SBATCH --output=cov_%j.out # Standard output and error log
#SBATCH --account=juannanzhou --qos=juannanzhou 
cd /orange/baer/briggsae
### step 5: Add read group information using picard:
module load samtools
out="individual_coverage_summary.txt"
echo -e "Sample\tSites_3x\tSites_10x" > "$out"

# If your BAMs are in the current directory:
for bam in *.markduplicates.bam; do
  s=$(basename "$bam" .bam)

  # count sites with depth >=3
  c3=$(samtools depth -a "$bam" | awk '$3>=3{n++} END{print n+0}')

  # count sites with depth >=10
  c10=$(samtools depth -a "$bam" | awk '$3>=10{n++} END{print n+0}')

  echo -e "${s}\t${c3}\t${c10}" >> "$out"
done

echo "Done – see $out"
