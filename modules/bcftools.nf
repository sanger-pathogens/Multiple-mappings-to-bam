process BCFTOOLS_CALL {
    tag "${meta.ID}"
    
    label "cpu_1"
    label "mem_500M"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/bcftools:1.11-c1'

    publishDir "${params.outdir}/${meta.ID}_${params.program}", mode: 'copy', overwrite: true, pattern: "*.{bcf,ploidy}"

    input:
    tuple val(meta), path(bam), path(mpileup)

    output:
    tuple val(meta), path(bam), path("${meta.ID}.bcf"), emit: called_ch
    path("${meta.ID}.ploidy")
    path("${meta.ID}_variant.bcf")
    path("${meta.ID}.bcf.csi")
    path("${meta.ID}_variant.bcf.csi")

    script:
    """
    echo "${meta.ID}    1" > ${meta.ID}.ploidy

    bcftools call -P ${params.prior} -O b -A -M -S ${meta.ID}.ploidy -${params.call} ${mpileup} > ${meta.ID}.bcf
    bcftools index ${meta.ID}.bcf
    bcftools call -P ${params.prior} -O b -A -M -v -S ${meta.ID}.ploidy -${params.call} ${mpileup} > ${meta.ID}_variant.bcf
    bcftools index ${meta.ID}_variant.bcf
    """
}