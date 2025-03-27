process INDEX_REF {
    tag "${ref}"

    label "cpu_1"
    label "mem_1"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    path(ref)

    output:
    tuple path(ref), path("${ref}.fai")

    script:
    """
    samtools faidx ${ref}
    """

}

process SAMTOOLS_SORT_BAM_AND_MAKE_HEADER {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.ID}.bam"), path("${meta.ID}.bam.bai"), emit: bam_ch
    tuple val(meta), path ("header.sam"), emit: header_ch

    script:
    """
    # Sort BAM file
    samtools sort -T ${meta.ID}.tmp ${bam} -o ${meta.ID}.bam

    # Index sorted BAM file
    samtools index ${meta.ID}.bam

    # Add read groups and fix header
    # sed is used to substitute the SO:unknown to SO:coordinate in the header
    samtools view -H ${meta.ID}.bam | sed 's/SO:unknown/SO:coordinate/g' | sed 's/\\\\x00//g' > header.sam
    """
}

process SAMTOOLS_SORT {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_1"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.ID}.bam")

    script:
    """
    samtools sort -T ${meta.ID}.tmp ${bam} -o ${meta.ID}.bam
    """
}

process SAMTOOLS_INDEX {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path (bam)

    output:
    tuple val(meta), path (bam), path("${meta.ID}.bam.bai"), emit: index_ch

    script:
    """
    samtools index ${bam}
    """

}

process SAMTOOLS_MERGE {
    tag "${meta.ID}"
    
    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    
    input:
    tuple val(meta), path(bam), path(bam_bai), path(header)

    output:
    tuple val(meta), path ("${meta.ID}.bam"), emit: bam_ch

    script:
    """
    samtools view -b -o tmphead.bam -H ${bam}
    samtools merge -c -p -f -r -h ${header} tmp.bam ${bam} tmphead.bam
    samtools reheader ${header} tmp.bam > ${meta.ID}.bam
    """
}

process SAMTOOLS_PILEUP {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    tuple val(meta), path(bam), path(bam_bai), path(ref), path(bwa_indexes)

    output:
    tuple val(meta), path(bam), path("${meta.ID}.mpileup")

    script:
    def dontuseanomolous = params.dontuseanomolous ? '' : ' -A '
    def BAQ = (params.BAQ == false) ? '' : ' -B '
    def overlaps = params.detectOverlaps ? params.detectOverlaps : '-x'

    def samtools_opts = "-t DP,DP4 -C 50 -L 1000 -d 1000 -m ${params.depth} ${dontuseanomolous} ${BAQ} ${overlaps} -ugf ${ref} ${bam}"
    """
    samtools mpileup $samtools_opts > ${meta.ID}.mpileup
    """

}
