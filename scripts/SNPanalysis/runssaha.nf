process runSsaha {
    input:
        val ref
        val fastqdir
        val name
        val pairedend
        val runname
        val ssahaquality
        val maxinsertsize
        val mininsertsize
        val circular
    
    output:

    script:
    if (!pairedend) {
        """
        ssaha2 -score ${ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -output sam_soft -outfile ${runname}/tmp1.sam ${ref} ${fastqdir}${name}.fastq"
        """
    } else {
        """
        ssaha2 -score ${ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -outfile ${runname}/tmp1.sam -pair ${mininsertsize},${maxinsertsize} -output sam_soft ${ref} ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq
        """
    }

    """
    samtools view -b -S ${runname}/tmp1.sam -t ${ref}.fai > ${runname}/tmp1.bam
    """

    if (pairedend && circular) {
        """
        fix_circular_bams.py -b ${runname}/tmp1.bam -o ${runname}/tmp
        rm ${runname}/tmp1.bam
        """
    } else {
        """
        mv ${runname}/tmp1.bam ${runname}/tmp.bam
        """
    }

    """
    samtools view -H ${runname}/tmp.bam > ${runname}/tmp2.sam
    cat ${runname}/tmp2.sam ${runname}/tmp1.sam > ${runname}/tmp.sam
    rm ${runname}/tmp2.sam ${runname}/tmp1.sam
    """
}