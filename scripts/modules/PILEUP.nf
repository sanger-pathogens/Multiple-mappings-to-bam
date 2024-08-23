process PILEUP {
    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    input:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai), path(ref)

    output:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai), path("${runname}/tmp.mpileup")

    script:
    runname = pools.runname
    """
    mkdir -p "${runname}"
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
        samtools mpileup -t DP,DP4 -C 50 -L 1000 -d 1000 -m ${params.depth} \$anomolous \$BAQ \$overlaps -ugf ${ref} ${name_bam} > ${runname}/tmp.mpileup
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
        samtools mpileup -t DP,DP4 -L 1000 -d 1000 -m ${params.depth} \$anomolous \$BAQ \$overlaps -ugf ${ref} ${name_bam} > ${runname}/tmp.mpileup
    fi
    """
}