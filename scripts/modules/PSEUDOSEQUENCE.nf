process PSEUDOSEQUENCE {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai), path(tmp_mpileup), path (name_bcf)

    output:
    path "${runname}/${name}*", emit: pseudosequence

    script:
    """
    mkdir -p "${runname}"
    if [ "${params.call}" = "m" ]; then
        bcf_2_pseudosequence.py -A -b ${name_bcf} -B ${name_bam} -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}
    elif [ "${params.call}" = "c" ]; then
        bcf_2_pseudosequence.py -A -b ${name_bcf} -B ${name_bam} -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}
    fi
    """
}