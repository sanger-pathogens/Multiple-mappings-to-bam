process SAMTOOLS_SORT {
    container 'quay.io/biocontainers/samtools:1.3--1'

    input:
    path tmp1_bam
    val name
    val runname

    output:
    path "${runname}/${name}..bam", emit: name_bam
    path "${runname}/tmphead.sam", emit: tmphead_sam
    path "${runname}/${name}.bam.bai", emit: bam_bai

    script:
    """
    mkdir -p ${runname}

    # Sort BAM file
    samtools sort ${tmp1_bam} > ${runname}/${name}.bam

    # Index sorted BAM file
    samtools index ${runname}/${name}.bam

    # Remove temporary BAM file
    rm ${tmp1_bam}

    # Add read groups and fix header
    # sed is used to substitute the SO:unknown to SO:coordinate in the header
    samtools view -H ${runname}/${name}.bam | sed 's/SO:unknown/SO:coordinate/g' | sed 's/\\\\x00//g' > ${runname}/tmphead.sam
    """
}

process SAMTOOLS_SORT_1 {
    container 'quay.io/biocontainers/samtools:1.3--1'

    input:
    path tmp1_bam

    output:
    path "${runname}/tmp.*", emit: tmp_bam

    script:
    """
    mkdir -p "${runname}"
    samtools sort ${tmp1_bam} > ${runname}/tmp.bam
    rm ${runname}/tmp1.bam
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
