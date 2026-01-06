process FORMAT_SSAHA_HEADER {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/void:0.0.1'

    stageInMode = 'copy'

    input:
    tuple val(meta), path(header)

    output:
    tuple val(meta), path(header)

    script:
    """
    now=\$(date +'%Y-%m-%dT%H:%M:%S')
    echo "@RG\tID:${meta.ID}\tCN:Sanger\tDT:"\$now"\tPG:SSAHA\tPL:ILLUMINA\tSM:${meta.ID}" >> ${header}
    """
}

process RUN_SSAHA {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/sangerpathogens/ssaha2:v2.5.5_cv3'

    input:
    tuple val(meta), path(name_1_fastq), path(name_2_fastq)
    tuple path(ref), path(ref_fai)

    output:
    tuple val(meta),  path("${final_name}"), path(ref_fai), emit: mapped_reads

    script:
    final_name = "${meta.ID}_mapped.sam"

    """
    ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -outfile ${final_name} -pair ${params.mininsertsize},${params.maxinsertsize} -output sam_soft ${ref} ${name_1_fastq} ${name_2_fastq}
    """
}

process FIX_CIRCULAR_BAMS {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/bam_filter:0.0.3'

    input:
    tuple val(meta),  path(bam)

    output:
    tuple val(meta),  path("${meta.ID}_fixed.bam"), emit: mapped_reads

    script:
    """
    fix_circular_bams.py -b ${bam} -o ${meta.ID}_fixed
    """
}