nextflow.enable.dsl=2

process RUN_SSAHA {
    container 'quay.io/sangerpathogens/ssaha2:v2.5.5_cv3'
    tag "${name}"
    input:
        val fastqdir
        val name
        val pairedend
        val runname

    script:
    """
        if [ "${pairedend}" = "false" ]; then
            ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -output sam_soft -outfile ${runname}/tmp1.sam ${ref} ${fastqdir}${name}.fastq
        else
            ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -outfile ${runname}/tmp1.sam -pair ${params.mininsertsize},${params.maxinsertsize} -output sam_soft ${ref} ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq
        fi

        samtools view -b -S ${runname}/tmp1.sam -t ${ref}.fai > ${runname}/tmp1.bam

        if [ "${pairedend}" = "true" ] && [ "${params.circular}" = "true" ]; then
            fix_circular_bams.py -b ${runname}/tmp1.bam -o ${runname}/tmp
            rm ${runname}/tmp1.bam
        else
            mv ${runname}/tmp1.bam ${runname}/tmp.bam
        fi

        samtools view -H ${runname}/tmp.bam > ${runname}/tmp2.sam
        cat ${runname}/tmp2.sam ${runname}/tmp1.sam > ${runname}/tmp.sam
        rm ${runname}/tmp2.sam ${runname}/tmp1.sam

    """
}
