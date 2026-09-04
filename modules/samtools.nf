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

process SAMTOOLS_SORT {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_1"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.ID}.bam"), emit: bam_ch

    script:
    """
    samtools sort -T ${meta.ID}.tmp ${bam} -o ${meta.ID}.bam
    """
}

process SAMTOOLS_INDEX {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_500M"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${meta.ID}.bam.bai"), emit: index_ch

    script:
    """
    samtools index ${bam}
    """

}

process UPDATE_HEADER_AND_RG_TAG {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_100M"
    label "time_30m"

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    tuple val(meta), path(sorted_bam)
    val(program)

    output:
    tuple val(meta), path(updated_bam), emit: updated_bam

    script:
    updated_bam = "${meta.ID}.bam"
    header = "header.sam"
    """
    # change SO:unknown to SO:coordinate in the header to indicate the BAM is coordinate sorted,
    # and strip stray null bytes
    samtools view -H ${sorted_bam} \\
    | sed -e 's/SO:unknown/SO:coordinate/g' -e 's/\\\\x00//g' > ${header}

    # GATK requires a readgroup (@RG) tag for indel realignment, so add one if missing
    if ! grep -q "^@RG" "${header}"; then
        now=\$(date +'%Y-%m-%dT%H:%M:%S')
        echo '@RG\tID:${meta.ID}\tCN:Sanger\tDT:'\$now'\tPG:${program}\tPL:ILLUMINA\tSM:${meta.ID}' >> ${header}
    fi

    samtools addreplacerg -r "\$(grep '^@RG' ${header})" -m overwrite_all ${sorted_bam} \\
    | samtools view -b - \\
    | samtools reheader ${header} - > ${updated_bam}
    """
}

process SAMTOOLS_SAM_TO_BAM {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_1"
    label "time_1"

    publishDir "${params.outdir}/${meta.ID}_${params.program}/raw_bams", mode: 'copy', overwrite: true, enabled: params.publish_raw_bam

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    tuple val(meta), path(sam), path(ref_fai)

    output:
    tuple val(meta), path("${meta.ID}.bam")

    script:
    """
    samtools view -b -S ${sam} -t ${ref_fai} > ${meta.ID}.bam
    """

}

process SAMTOOLS_PILEUP {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_250M"
    label "time_1"
    
    publishDir "${params.outdir}/${meta.ID}_${params.program}", mode: 'copy', overwrite: true, pattern: "*.mpileup"

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    tuple val(meta), path(bam), path(bam_bai), path(ref)

    output:
    tuple val(meta), path(bam), path("${meta.ID}.mpileup")

    script:
    def dontuseanomolous = params.dontuseanomolous ? '' : ' -A '
    def BAQ = (params.BAQ == false) ? '' : ' -B '
    def overlaps = params.detectOverlaps ? '' : '-x'

    def samtools_opts = "-t DP,DP4 -C 50 -L 1000 -d 1000 -m ${params.depth} ${dontuseanomolous} ${BAQ} ${overlaps} -ugf ${ref} ${bam}"
    """
    samtools mpileup $samtools_opts > ${meta.ID}.mpileup
    """

}
