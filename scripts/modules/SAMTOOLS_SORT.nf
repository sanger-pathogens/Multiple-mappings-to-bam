process SAMTOOLS_SORT {
    container 'samtools-1.3'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pools), path(file1), path(file2), path(tmp1_bam), val(cmdline)

    output:
    tuple val(pools), path(file1), path(file2), path ("${runname}/${name}.bam"), val(cmdline), path ("${runname}/tmphead.sam"), path("${runname}/${name}.bam.bai")

    script:
    runname = pools.runname
    name = pools.name
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
    container 'samtools-1.3'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pools), path(file1), path(file2), path (tmp1_bam), path(tmphead_sam), path(bam_bai)

    output:
    tuple val(pools), path(file1), path(file2), path ("${runname}/tmp.bam"), path(tmphead_sam), path(bam_bai)

    script:
    """
    mkdir -p "${runname}"
    samtools sort ${tmp1_bam} > ${runname}/tmp.bam
    rm ${runname}/tmp1.bam
    """
}

process SAMTOOLS_INDEX {
    container 'samtools-1.3'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai)

    output:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path("*.bai")

    script:
    """
    samtools index ${name_bam}
    """

}
process SAMTOOLS_MERGE {
    container 'samtools-1.3'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    
    input:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai)

    output:
    tuple val(pools), path(file1), path(file2), path ("${runname}/tmp1.bam"), path(tmphead_sam), path(bam_bai)
    path ("${runname}/tmphead.bam")

    script:
    runname = pools.runname
    name = pools.name
    """
    mkdir -p ${runname}
    samtools view -b -o ${runname}/tmphead.bam -H ${name_bam}
    samtools merge -c -p -f -r -h ${tmphead_sam} ${runname}/tmp.bam ${name_bam} ${runname}/tmphead.bam
    mv ${runname}/tmp.bam ${runname}/tmp1.bam
    rm ${name}.bam
    """
}
