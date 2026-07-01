module load bcftools

module load htslib



for gvcf in *.dv.g.vcf.gz; do

    sample="${gvcf%.dv.g.vcf.gz}"          # e.g. HK_204.hifi_reads

    echo "$sample" > tmp_sample.txt

    bcftools reheader -s tmp_sample.txt "$gvcf" -o "${gvcf%.vcf.gz}.renamed.vcf.gz"

    bcftools index -t "${gvcf%.vcf.gz}.renamed.vcf.gz"

    echo "[DONE] $gvcf → sample name: $sample"

done

rm tmp_sample.txt



# Build new list pointing to renamed files

ls *.dv.g.renamed.vcf.gz > gvcfs.fixed.list
