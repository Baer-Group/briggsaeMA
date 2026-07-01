for bam in *.hifi_reads.bam; do

    echo "===== $bam =====" >> mapping_summary.txt

    samtools flagstat "$bam" >> mapping_summary.txt

    echo "" >> mapping_summary.txt

done

for bam in *.hifi_reads.bam; do

    sample=${bam%.hifi_reads.bam}

    samtools flagstat "$bam" | awk -v s="$sample" '/mapped \(/ {print s,$5}'

done > mapping_rate.txt
