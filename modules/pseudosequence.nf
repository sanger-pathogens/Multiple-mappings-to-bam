process PSEUDOSEQUENCE {
    tag "${meta.ID}"
    
    label "cpu_1"
    label "mem_500M"
    label "time_1"
    
    publishDir "${params.outdir}/${meta.ID}_${params.program}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/bcf_2_pseudosequence:0.0.2'
    
    input:
    tuple val(meta), path(bam), path(bcf)

    output:
    tuple val(meta), path("${meta.ID}.mfa"), emit: pseudosequence

    script:
    """
    if [ "${params.call}" = "m" ]; then
        bcf_2_pseudosequence.py -A -b ${bcf} -B ${bam} -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${meta.ID}
    elif [ "${params.call}" = "c" ]; then
        bcf_2_pseudosequence.py -A -b ${bcf} -B ${bam} -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${meta.ID}
    fi
    """
}