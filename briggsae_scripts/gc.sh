#!/bin/bash

module load samtools

echo -e "Sample\tGC_percent" > GC_content_per_sample.tsv



for bam in *.markduplicates.bam

do

    sample=$(basename $bam .bam)



    gc=$(samtools stats $bam | awk '/^GC/ {print $3}')



    echo -e "${sample}\t${gc}" >> GC_content_per_sample.tsv

done
