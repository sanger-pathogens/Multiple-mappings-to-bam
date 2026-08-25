process FORMAT_BWA_HEADER {
    tag "${meta.ID}"
    
    label "cpu_1"
    label "mem_100M"
    label "time_30m"

    container 'quay.io/ssd28/gsoc-experimental/void:0.0.1'

    stageInMode = 'copy'

    input:
    tuple val(meta), path(header)

    output:
    tuple val(meta), path(header), emit: header_ch

    script:
    """
    now=\$(date +'%Y-%m-%dT%H:%M:%S')
    echo '@RG\tID:${meta.ID}\tCN:Sanger\tDT:'\$now'\tPG:BWA MEM\tPL:ILLUMINA\tSM:${meta.ID}' >> ${header}
    """
}


process RUN_BWA {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_1"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/run-bwa:0.0.2'

    publishDir "${params.outdir}/${meta.ID}_${params.program}/raw_bams", mode: 'copy', overwrite: true, saveAs: { filename -> "${meta.ID}.bam"}, enabled: params.publish_raw_bam

    input:
    tuple val(meta), path(fastq1), path(fastq2)
    tuple path(ref), path(bwa_index_files)
    
    output:
    tuple val(meta), path ("mapped.bam"), emit: mapped_ch

    script:
    """
    bwa mem -v 1 -M -a -t 1 ${ref} ${fastq1} ${fastq2} > mapped.sam
    samtools view -b -S mapped.sam -t *.fai > mapped.bam
    rm -f tmp.sam
    """
}

process BWA_INDEX {
    label "cpu_1"
    label "mem_250M"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/run-bwa:0.0.2'

    input:
    path(ref)

    output:
    tuple path(ref), path("${ref}.*"), emit: index_ch

    script:
    """
    bwa index ${ref}
    samtools faidx ${ref}
    """
}

