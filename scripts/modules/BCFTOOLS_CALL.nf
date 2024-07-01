process BCFTOOLS_CALL {
    input:
    val runname
    val name
    path tmp_mpileup

    output:
    path "${runname}/${name}.ploidy", emit: name_ploidy
    path "${runname}/${name}.bcf", emit: name_bcf
    path "${runname}/${name}_variant.bcf", emit: name_variant_bcf
    path "${runname}/${name}.bcf.csi", emit: name_bcf
    path "${runname}/${name}_variant.bcf", emit: name_variant_bcf_csi

    script:
    """
    mkdir -p ${runname}
    echo "${name}    1" > ${runname}/${name}.ploidy

    bcftools call -P ${params.prior} -O b -A -M -S ${runname}/${name}.ploidy -${params.call} ${tmp_mpileup} > ${runname}/${name}.bcf
    bcftools index ${runname}/${name}.bcf
    bcftools call -P ${params.prior} -O b -A -M -v -S ${runname}/${name}.ploidy -${params.call} ${tmp_mpileup} > ${runname}/${name}_variant.bcf
    bcftools index ${runname}/${name}_variant.bcf
    """
}