process SAMTOOLS_SORT {
    container 'quay.io/biocontainers/samtools:1.3--1'

    input:
    tuple val(name), path(tmp1_bam)

    output:
    tuple val(name), path ("${params.runname}/${name}.bam"), emit: name_bam
    tuple val(name), path ("${params.runname}/tmphead.sam"), emit: tmphead_sam
    path "${params.runname}/${name}.bam.bai", emit: bam_bai

    script:
    """
    mkdir -p ${params.runname}

    # Sort BAM file
    samtools sort ${tmp1_bam} > ${params.runname}/${name}.bam

    # Index sorted BAM file
    samtools index ${params.runname}/${name}.bam

    # Remove temporary BAM file
    rm ${tmp1_bam}

    # Add read groups and fix header
    # sed is used to substitute the SO:unknown to SO:coordinate in the header
    samtools view -H ${params.runname}/${name}.bam | sed 's/SO:unknown/SO:coordinate/g' | sed 's/\\\\x00//g' > ${params.runname}/tmphead.sam
    """
}

process SAMTOOLS_SORT_1 {
    container 'quay.io/biocontainers/samtools:1.3--1'

    input:
    tuple val(name), path (tmp1_bam)

    output:
    tuple val(name), path ("${params.runname}/tmp.bam"), emit: tmp_bam

    script:
    """
    mkdir -p "${params.runname}"
    samtools sort ${tmp1_bam} > ${params.runname}/tmp.bam
    rm ${params.runname}/tmp1.bam
    """
}

process SAMTOOLS_INDEX {
    container 'quay.io/biocontainers/samtools:1.3--1'

    input:
    path name_bam

    output:
    path "*.bai", emit: bam_bai

    script:
    """
    samtools index ${name_bam}
    """

}
process SAMTOOLS_MERGE {
    input:
    tuple val(name), path (name_bam)
    path tmphead_sam

    output:
    tuple val(name), path ("${params.runname}/tmp1.bam"), emit: tmp1_bam
    path "${params.runname}/tmphead.bam", emit; tmphead_bam

    script:
    """
    mkdir -p ${params.runname}
    samtools view -b -o ${params.runname}/tmphead.bam -H ${name_bam}
    samtools merge -c -p -f -r -h ${tmphead_sam} ${params.runname}/tmp.bam ${name_bam} ${params.runname}/tmphead.bam
    mv ${params.runname}/tmp.bam ${params.runname}/tmp1.bam
    rm ${params.runname}/${name}.bam
    """
}
