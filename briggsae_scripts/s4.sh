#!/bin/sh
module load samtools/

module load bamtools/

for i in *.sam

do

 samtools view -b -q 3 -o "${i/%.sam/.bam}" "$i"

done

rm *.sam

for i in *.bam

do

 samtools sort -O bam -T tmp_ -o "${i/%.bam/.output.bam}" "$i"

done

for i in *.output.bam

do

 bamtools filter -isProperPair true -in "$i" -out "${i/%.output.bam/.PP.bam}"

done

rm *.output.bam
