process FORMAT_SMALT_HEADER {
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
    echo "@RG\tID:${meta.ID}\tCN:Sanger\tDT:"\$now"\tPG:SMALT\tPL:ILLUMINA\tSM:${meta.ID}" >> ${header}
    """
}


process RUN_SMALT {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/run-smalt:0.0.2'
    
    input:
    tuple val(meta), path(name_1_fastq), path(name_2_fastq)
    tuple path(ref), path(ref_sma), path(ref_smi)

    output:
    tuple val(meta), path(final_name), emit: mapped_reads

    script:
    final_name = "${meta.ID}_mapped.bam"

    allow_multimapping = params.allow_multimapping ? '-r 0' : '-r -1'

    """
    smalt map -y ${params.nomapid} \\
        -x \\
        ${allow_multimapping} \\
        -i ${params.maxinsertsize} \\
        -j ${params.mininsertsize} \\
        -f bam \\
        -o ${final_name} \\
        ${ref_sma.baseName} \\
        ${name_1_fastq} \\
        ${name_2_fastq}
    """
}

process SMALT_INDEX {
    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/run-smalt:0.0.2'

    input:
    path(ref)

    output:
    tuple path(ref), path("${ref.baseName}_index.sma"), path("${ref.baseName}_index.smi")

    script:
    """
    if [ "${params.human}" == "True" ]; then
        smalt index -k 20 -s 13 ${ref.baseName}_index ${ref}
    else 
        smalt index -k 13 -s 1 ${ref.baseName}_index ${ref}
    fi
    """
}