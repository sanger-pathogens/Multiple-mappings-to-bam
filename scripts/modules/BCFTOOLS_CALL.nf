process BCFTOOLS_CALL {
    label "cpu_1"
    label "mem_16"
    label "time_1"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/bcftools:1.11'

    input:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai), path(tmp_mpileup)

    output:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai), path(tmp_mpileup), path ("${runname}/${name}.bcf")
    path "${runname}/${name}.ploidy", emit: name_ploidy
    path "${runname}/${name}_variant.bcf", emit: name_variant_bcf
    path "${runname}/${name}.bcf.csi", emit: name_bcf_csi
    path "${runname}/${name}_variant.bcf", emit: name_variant_bcf_csi

    script:
    runname = pools.runname
    name = pools.name
    """
    mkdir -p ${runname}
    echo "${name}    1" > ${runname}/${name}.ploidy

    bcftools call -P ${params.prior} -O b -A -M -S ${runname}/${name}.ploidy -${params.call} ${tmp_mpileup} > ${runname}/${name}.bcf
    bcftools index ${runname}/${name}.bcf
    bcftools call -P ${params.prior} -O b -A -M -v -S ${runname}/${name}.ploidy -${params.call} ${tmp_mpileup} > ${runname}/${name}_variant.bcf
    bcftools index ${runname}/${name}_variant.bcf
    """
}