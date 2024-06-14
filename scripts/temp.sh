smalt map -y 0 -x -r 0 -i 500 -j 100 -f bam -o run_name/tmp1.bam temp_name.index path/to/fastq/sample_name_1.fastq path/to/fastq/sample_name_2.fastq
samtools view -b -S run_name/tmp1.sam -t path/to/reference.fai > run_name/tmp1.bam
rm run_name/tmp1.sam
