process INDEX {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val tmpname
    path ref

    output:
    path "${ref}*", emit: mapping_index_ch

    script:
    """
    if [ "${params.program}" == "BWA" ]; then
        bwa index ${ref}
    elif [ "${params.program}" == "SMALT" ]; then
        if [ "${params.human}" == "True" ]; then
            smalt index -k 20 -s 13 ${tmpname}.index ${ref}
        else 
            smalt index -k 13 -s 1 ${tmpname}.index ${ref}
        fi
    fi

    samtools faidx ${ref}
    """
}