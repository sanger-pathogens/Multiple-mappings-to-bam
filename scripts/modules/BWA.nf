process BWA {
    input:
    tuple val(pools), path(file1), path(file2), path (name_bam), val(cmdline), path(tmphead_sam), path(bam_bai)

    output:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai)

    script:
    name = pools.name
    """
    echo '@RG\\tID:${name}\\tCN:Sanger\\tDT:\$now\\tPG:BWA MEM\\tPL:ILLUMINA\\tSM:${name}' >> ${tmphead_sam}
    """
}


process RUN_BWA {
    container 'bwa-samtools:latest'
    tag "${name}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pools), path(name_1_fastq), path(name_2_fastq), path(ref), path(ref_fai)
    tuple path(f1), path(f2), path(f3), path(f4), path(f5), path(f6)
    
    output:
    tuple val(pools), path(name_1_fastq), path(name_2_fastq), path ("${runname}/tmp1.bam"), val(cmdline)

    script:
    name = pools.name
    runname = pools.runname
    is_zipped = pools.is_zipped
    pairedend = pools.pairedend
    cmdline = ""
    
    """
    mkdir -p "${runname}"
    if [ "${is_zipped}" = "true" ]; then
        if [ "${pairedend}" = "true" ]; then
            bwa mem -v 1 -M -a -t 1 ${ref} ${name_1_fastq} ${name_2_fastq} > ${runname}/tmp.sam
        else
            bwa mem -v 1 -M -a -t 1 ${ref} ${name_1_fastq} > ${runname}/tmp.sam
        fi
    else
        if [ "${pairedend}" = "true" ]; then
            bwa mem -v 1 -M -a -t 1 ${ref} ${name_1_fastq} ${name_2_fastq} > ${runname}/tmp.sam
        else
            bwa mem -v 1 -M -a -t 1 ${ref} ${name_1_fastq} > ${runname}/tmp.sam
        fi
    fi

    samtools view -b -S ${runname}/tmp.sam -t ${ref_fai} > ${runname}/tmp1.bam
    rm -f ${runname}/tmp.sam
    """ 
}

process BWA_INDEX {
    container 'bwa-samtools:latest'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ref

    output:
    tuple path("*.amb"), path("*.ann"), path("*.bwt"), path("*.pac"), path("*.sa"), path("*.alt")
    path("*.fai"), emit: fai

    when:
    params.program == "BWA"

    script:
    """
    bwa index ${ref}
    samtools faidx ${ref}
    touch ${ref}.alt
    """
}

