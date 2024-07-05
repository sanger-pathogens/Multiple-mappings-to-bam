process PSEUDOSEQUENCE {
    input:
    path name_bcf
    path name_bam

    output:
    path "${params.runname}/${name}*", emit: pseudosequence

    script:
    """
    mkdir -p "${params.runname}"
    if [ "${params.call}" = "m" ]; then
        bcf_2_pseudosequence.py -A -b ${params.runname}/${name}.bcf -B ${params.runname}/${name}.bam -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${params.runname}/${name}
    elif [ "${params.call}" = "c" ]; then
        bcf_2_pseudosequence.py -A -b ${params.runname}/${name}.bcf -B ${params.runname}/${name}.bam -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${params.runname}/${name}
    fi
    """
}