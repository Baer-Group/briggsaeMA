#!/usr/bin/env bash

set -euo pipefail



# Pair L001 vs L002 markduplicates BAMs

# Example:

#   202_POOLRET91_S9_L001.markduplicates.bam  <->  202_POOLRET91_S9_L002.markduplicates.bam



out="pairs.txt"

: > "$out"  # truncate/create



shopt -s nullglob



# Loop over all L001 BAMs and find matching L002 BAM by filename stem

for f1 in *_L001.markduplicates.bam; do

  f2="${f1/_L001./_L002.}"



  if [[ -f "$f2" ]]; then

    printf "%s\t%s\n" "$f1" "$f2" >> "$out"

  else

    # If you prefer to skip unmatched, replace this block with: continue

    printf "%s\t\n" "$f1" >> "$out"

  fi

done



echo "Done — wrote $(wc -l < "$out") rows to $out"


