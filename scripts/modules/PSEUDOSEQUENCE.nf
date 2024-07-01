process PSEUDOSEQUENCE {
    input:
    path name_bcf
    path name_bam

    output:
    path "${runname}/${name}*", emit: pseudosequence

    script:
    """
    mkdir -p "${runname}"
    if [ "${params.call}" = "m" ]; then
        bcf_2_pseudosequence.py -A -b ${runname}/${name}.bcf -B ${runname}/${name}.bam -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}
    elif [ "${params.call}" = "c" ]; then
        bcf_2_pseudosequence.py -A -b ${runname}/${name}.bcf -B ${runname}/${name}.bam -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}
    fi
    """
}