nextflow.enable.dsl=2

process runSsaha {
    input:
        path bashfile
        val fastqdir
        val name
        val pairedend
        val runname
    
    output:
        path bashfile

    script:
    """
        echo "\nRunning Ssaha on ${name}..."
        touch ${bashfile}

        if [ "${pairedend}" = "false" ]; then
            echo "ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -output sam_soft -outfile ${runname}/tmp1.sam ${ref} ${fastqdir}${name}.fastq" >> ${bashfile}
        else
            echo "ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -outfile ${runname}/tmp1.sam -pair ${params.mininsertsize},${params.maxinsertsize} -output sam_soft ${ref} ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq" >> ${bashfile}
        fi

        echo "samtools view -b -S ${runname}/tmp1.sam -t ${ref}.fai > ${runname}/tmp1.bam" >> ${bashfile}

        if [ "${pairedend}" = "true" ] && [ "${params.circular}" = "true" ]; then
            echo "fix_circular_bams.py -b ${runname}/tmp1.bam -o ${runname}/tmp" >> ${bashfile}
            echo "rm ${runname}/tmp1.bam" >> ${bashfile}
        else
            echo "mv ${runname}/tmp1.bam ${runname}/tmp.bam" >> ${bashfile}
        fi

        echo "samtools view -H ${runname}/tmp.bam > ${runname}/tmp2.sam" >> ${bashfile}
        echo "cat ${runname}/tmp2.sam ${runname}/tmp1.sam > ${runname}/tmp.sam" >> ${bashfile}
        echo "rm ${runname}/tmp2.sam ${runname}/tmp1.sam" >> ${bashfile}

    """
}