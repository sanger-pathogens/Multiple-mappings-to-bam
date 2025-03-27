process MARK_DUPLICATES {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/picard:1.126'
    
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.ID}_marked.bam"), emit: deduped_ch
    path("${meta.ID}_metrics.txt")

    script:
    """
    java -jar /usr/local/bin/picard.jar MarkDuplicates \\
        INPUT=${bam_file} \\
        OUTPUT=${meta.ID}_marked.bam \\
        METRICS_FILE=${meta.ID}_metrics.txt
    """
}

process SEQUENCE_DICT {
    tag "${ref}"

    label "cpu_1"
    label "mem_500M"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/picard:1.126'
    
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.txt"

    input:
    tuple path(ref), path(fai)

    output:
    tuple path(ref), path(fai) , path(dict), emit: ref_ch

    script:
    dict="${ref.simpleName}.dict"
    """
    java -jar /usr/local/bin/picard.jar CreateSequenceDictionary \\
        R=${ref} \\
        O=${dict}
    """
}
