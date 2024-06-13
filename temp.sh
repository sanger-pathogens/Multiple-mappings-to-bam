ssaha2 -score 30 -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -output sam_soft -outfile /tmp1.sam /home/ref.txt name.fastq
samtools view -b -S /tmp1.sam -t /home/ref.txt.fai > /tmp1.bam
fix_circular_bams.py -b /tmp1.bam -o /tmp
rm /tmp1.bam
samtools view -H /tmp.bam > /tmp2.sam
cat /tmp2.sam /tmp1.sam > /tmp.sam
rm /tmp2.sam /tmp1.sam
