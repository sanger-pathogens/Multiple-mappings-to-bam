process RUN_SSAHA {
    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/sangerpathogens/ssaha2:v2.5.5_cv3'
    tag "${name}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(pools), path(name_1_fastq), path(name_2_fastq), path(ref), path(ref_fai)

    output:
    tuple val(pools), path(name_1_fastq), path(name_2_fastq), path("${runname}/tmp1.bam"), val(cmdline)

    script:
    runname = pools.runname
    pairedend = pools.pairedend
    cmdline=""
    """
    mkdir -p "${runname}"
    if [ "${pairedend}" = "false" ]; then
        ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -output sam_soft -outfile ${runname}/tmp1.sam ${ref} ${name_fastq}
    else
        ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -outfile ${runname}/tmp1.sam -pair ${params.mininsertsize},${params.maxinsertsize} -output sam_soft ${ref} ${name_1_fastq} ${name_2_fastq}
    fi

    samtools view -b -S ${runname}/tmp1.sam -t ${ref_fai} > ${runname}/tmp1.bam

    if [ "${pairedend}" = "true" ] && [ "${params.circular}" = "true" ]; then
        fix_circular_bams.py -b ${runname}/tmp1.bam -o ${runname}/tmp
        rm ${runname}/tmp1.bam
    else
        mv ${runname}/tmp1.bam ${runname}/tmp.bam
    fi

    samtools view -H ${runname}/tmp.bam > ${runname}/tmp2.sam
    cat ${runname}/tmp2.sam ${runname}/tmp1.sam > ${runname}/tmp.sam
    samtools view -b -S ${runname}/tmp.sam -t ${ref_fai} > ${runname}/tmp1.bam
    rm -f ${runname}/tmp.sam
    rm ${runname}/tmp2.sam ${runname}/tmp1.sam
    """
}
