process RUN_SSAHA {
    container 'quay.io/sangerpathogens/ssaha2:v2.5.5_cv3'
    tag "${name}"

    input:
    path ref
    tuple val(name), path (name_fastq)  //optional
    path name_1_fastq                   //optional
    path name_2_fastq                   //optional
    path ref_fai
    val fastqdir
    val pairedend

    output:
    path "${params.runname}/tmp*", emit: tmp_bam
    path "${params.runname}/tmp.sam", emit: tmp_sam

    script:
    """
    mkdir -p "${runname}"
    if [ "${pairedend}" = "false" ]; then
        ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -output sam_soft -outfile ${params.runname}/tmp1.sam ${ref} ${name_fastq}
    else
        ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -outfile ${params.runname}/tmp1.sam -pair ${params.mininsertsize},${params.maxinsertsize} -output sam_soft ${ref} ${name_1_fastq} ${name_2_fastq}
    fi

    samtools view -b -S ${params.runname}/tmp1.sam -t ${ref_fai} > ${params.runname}/tmp1.bam

    if [ "${pairedend}" = "true" ] && [ "${params.circular}" = "true" ]; then
        fix_circular_bams.py -b ${params.runname}/tmp1.bam -o ${params.runname}/tmp
        rm ${params.runname}/tmp1.bam
    else
        mv ${params.runname}/tmp1.bam ${params.runname}/tmp.bam
    fi

    samtools view -H ${params.runname}/tmp.bam > ${params.runname}/tmp2.sam
    cat ${params.runname}/tmp2.sam ${params.runname}/tmp1.sam > ${params.runname}/tmp.sam
    rm ${params.runname}/tmp2.sam ${params.runname}/tmp1.sam
    """
}
