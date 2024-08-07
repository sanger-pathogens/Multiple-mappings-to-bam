process PSEUDOSEQUENCE {
    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/bcf_2_pseudosequence:0.0.1'
    
    input:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai), path(tmp_mpileup), path (name_bcf)

    output:
    path "${runname}/${name}*", emit: pseudosequence

    script:
    runname = pools.runname
    name = pools.name
    """
    mkdir -p "${runname}"
    if [ "${params.call}" = "m" ]; then
        python2 /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py -A -b ${name_bcf} -B ${name_bam} -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}
    elif [ "${params.call}" = "c" ]; then
        python2 /opt/bcf_2_pseudosequence/bcf_2_pseudosequence.py -A -b ${name_bcf} -B ${name_bam} -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}
    fi
    """
}