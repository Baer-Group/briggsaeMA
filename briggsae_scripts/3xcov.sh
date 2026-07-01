#!/bin/sh

#SBATCH --job-name=3xcov # Job name

#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)

#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail

#SBATCH --cpus-per-task=4 # Number of CPU cores per task

#SBATCH --ntasks=1 # Run a single task

#SBATCH --mem=20gb # Memory limit

#SBATCH --time=80:00:00 # Time limit hrs:min:sec

#SBATCH --output=3xcov_%j.out # Standard output and error log

#SBATCH --account=baer --qos=baer


# Load the necessary module for samtools

module load samtools


# Output file

out="common_3x_summary.txt"

# Header

echo -e "Sample1\tSample2\tCommonSites_3x" > "$out"



# Loop over each pair (tab‐separated: Efile<TAB>Gfile)

while IFS=$'\t' read -r bam1 bam2; do

  # Skip pairs with no partner

  if [[ -z "$bam2" ]]; then

    echo -e "$(basename "$bam1")\tNA\t0" >> "$out"

    continue

  fi



  # Basenames (no .bam)

  s1=$(basename "$bam1" .bam)

  s2=$(basename "$bam2" .bam)



  # Count positions where both depths ≥ 3

  shared=$(samtools depth -a "$bam1" "$bam2" | awk '$3>=3 && $4>=3' | wc -l)



  # Append to summary

  echo -e "${s1}\t${s2}\t${shared}" >> "$out"

done < pairs.txt



echo "Done – see $out"


