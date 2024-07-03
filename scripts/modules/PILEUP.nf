process PILEUP {
    input:
    path ref
    tuple val(name), path (name_bam)

    output:
    tuple val(name), path ("${params.runname}/tmp.mpileup"), emit: tmp_mpileup

    script:
    """
    mkdir -p "${params.runname}"
    if [ ${params.anomolous} ]
    then
        anomolous=" -A "
    else
        anomolous=""
    fi

    if [ "${params.program}" = "bwa" ] || [ "${params.program}" = "BWA" ]
    then
        if [ ${params.BAQ} = "false" ]
        then
            BAQ="-B"
        else
            BAQ=""
        fi
        if [ ${params.detectOverlaps} ]
        then
            overlaps=""
        else
            overlaps="-x"
        fi
        samtools mpileup -t DP,DP4 -C 50 -L 1000 -d 1000 -m ${params.depth} \$anomolous \$BAQ \$overlaps -ugf ${ref} ${name_bam} > ${params.runname}/tmp.mpileup
    else
        if [ ${params.BAQ} = "true" ]
        then
            BAQ=""
        else
            BAQ="-B"
        fi
        if [ ${params.detectOverlaps} = "true" ]
        then
            overlaps="-x"
        else
            overlaps=""
        fi
        samtools mpileup -t DP,DP4 -L 1000 -d 1000 -m ${params.depth} \$anomolous \$BAQ \$overlaps -ugf ${ref} ${name_bam} > ${params.runname}/tmp.mpileup
    fi
    """
}