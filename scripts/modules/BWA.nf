process BWA {
    input:
    tuple val(name), path(tmphead_sam)

    output:
    tuple val(name), path(tmphead_sam)

    script:
    """
    echo '@RG\\tID:${name}\\tCN:Sanger\\tDT:\$now\\tPG:BWA MEM\\tPL:ILLUMINA\\tSM:${name}' >> ${tmphead_sam}
    """
}


process RUN_BWA {
    container 'quay.io/ssd28/gsoc-experimental/run-bwa:0.0.1'
    tag "${name}"

    input:
    path ref
    tuple val(name), path name_fastq
    path name_1_fastq
    path name_2_fastq
    path ref_fai
    val fastqdir
    val pairedend
    val is_zipped
    
    output:
    path "${params.runname}/tmp1.bam", emit: tmp1_bam

    script:
    """
    mkdir -p "${runname}"
    if [ "${is_zipped}" = "true" ]; then
        if [ "${pairedend}" = "true" ]; then
            bwa mem -v 1 -M -a -t 1 ${ref} ${name_1_fastq} ${name_2_fastq} > ${params.runname}/tmp.sam
        else
            bwa mem -v 1 -M -a -t 1 ${ref} ${name_fastq} > ${params.runname}/tmp.sam
        fi
    else
        if [ "${pairedend}" = "true" ]; then
            bwa mem -v 1 -M -a -t 1 ${ref} ${name_1_fastq} ${name_2_fastq} > ${params.runname}/tmp.sam
        else
            bwa mem -v 1 -M -a -t 1 ${ref} ${name_fastq} > ${params.runname}/tmp.sam
        fi
    fi

    samtools view -b -S ${params.runname}/tmp.sam -t ${ref_fai} > ${params.runname}/tmp1.bam
    rm -f ${params.runname}/tmp.sam
    """
}
