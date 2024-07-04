process BCFTOOLS_CALL {
    input:
    tuple val(name), path (tmp_mpileup)

    output:
    path "${params.runname}/${name}.ploidy", emit: name_ploidy
    path "${params.runname}/${name}.bcf", emit: name_bcf
    path "${params.runname}/${name}_variant.bcf", emit: name_variant_bcf
    path "${params.runname}/${name}.bcf.csi", emit: name_bcf_csi
    path "${params.runname}/${name}_variant.bcf", emit: name_variant_bcf_csi

    script:
    """
    mkdir -p ${params.runname}
    echo "${name}    1" > ${params.runname}/${name}.ploidy

    bcftools call -P ${params.prior} -O b -A -M -S ${params.runname}/${name}.ploidy -${params.call} ${tmp_mpileup} > ${params.runname}/${name}.bcf
    bcftools index ${params.runname}/${name}.bcf
    bcftools call -P ${params.prior} -O b -A -M -v -S ${params.runname}/${name}.ploidy -${params.call} ${tmp_mpileup} > ${params.runname}/${name}_variant.bcf
    bcftools index ${params.runname}/${name}_variant.bcf
    """
}