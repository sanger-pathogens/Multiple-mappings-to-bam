process FORMAT_SMALT_HEADER {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/void:0.0.1'

    stageInMode = 'copy'

    input:
    tuple val(meta), path(header)

    output:
    tuple val(meta), path(header)

    script:
    newsmalt = "true"
    """
    now=\$(date +'%Y-%m-%dT%H:%M:%S')
    echo "@RG\tID:${meta.ID}\tCN:Sanger\tDT:"\$now"\tPG:SMALT\tPL:ILLUMINA\tSM:${meta.ID}" >> ${tmphead_sam}
    if [ ${params.domapping} ] && [ ${newsmalt} = "false" ]
    then
        smaltversion=\$( smalt version | grep Version | awk '{print \$2}' )
        echo "@PG\tID:SMALT\tPN:SMALT\tCL:${cmdline}\tVN:\$smaltversion" >> ${tmphead_sam}
    fi
    """
}


process RUN_SMALT {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/run-smalt:0.0.2'
    
    input:
    tuple val(meta), path(name_1_fastq), path(name_2_fastq)
    tuple path(ref), path(smalt_indexes)

    output:
    tuple val(meta), path(name_1_fastq), path(name_2_fastq), path("tmp1.bam"), env(cmdline)


    script:
    newsmalt = true
    domapping = meta.domapping
    pairedend = meta.pairedend
    fastqdir = meta.fastqdir
    bam = meta.bam
    """
    smaltoutput="bam"
    smaltoutputsuffix="bam"
    rbit=""
    cmdline=""

    if [ "${domapping}" = "true" ]; then

        if [ "${newsmalt}" = "true" ]; then
            smaltoutput="bam"
            smaltoutputsuffix="bam"
        else
            smaltoutput="samsoft"
            smaltoutputsuffix="sam"
        fi

        if [ "${pairedend}" = "true" ]; then
            if [ "${params.maprepeats}" = "true" ]; then
                smalt map -y ${params.nomapid} -x -r 0 -i ${params.maxinsertsize} -j ${params.mininsertsize} -f \$smaltoutput -o tmp1.\$smaltoutputsuffix ${meta.ID}_tmp.index ${name_1_fastq} ${name_2_fastq}
                cmdline="map -y ${params.nomapid} -x -r 0 -i ${params.maxinsertsize} -j ${params.mininsertsize} -f \$smaltoutput -o tmp1.\$smaltoutputsuffix ${meta.ID}_tmp.index ${fastqdir}${meta.ID}_1.fastq ${fastqdir}${meta.ID}_2.fastq"
            else
                if [ "${newsmalt}" = "true" ]; then
                    rbit=" -r -1"
                else
                    rbit=""
                fi
                smalt map -y ${params.nomapid}\$rbit -x -i ${params.maxinsertsize} -j ${params.mininsertsize} -f \$smaltoutput -o tmp1.\$smaltoutputsuffix ${meta.ID}_tmp.index ${name_1_fastq} ${name_2_fastq}
                cmdline="map -y ${params.nomapid}\$rbit -x -i ${params.maxinsertsize} -j ${params.mininsertsize} -f \$smaltoutput -o tmp1.\$smaltoutputsuffix ${meta.ID}_tmp.index ${fastqdir}${meta.ID}_1.fastq ${fastqdir}${meta.ID}_2.fastq"
            fi
        else
            if [ "${params.maprepeats}" = "true" ]; then
                smalt map -y ${params.nomapid} -x -r 0 -f \$smaltoutput -o tmp1.\$smaltoutputsuffix ${meta.ID}_tmp.index ${name_1_fastq}
                cmdline="map -y ${params.nomapid} -x -r 0 -f \$smaltoutput -o tmp1.\$smaltoutputsuffix ${meta.ID}_tmp.index ${fastqdir}${meta.ID}.fastq"
            else
                if [ "${newsmalt}" = "true" ]; then
                    \$rbit=" -r -1"
                else
                    \$rbit=""
                fi
                smalt map -y ${params.nomapid}\$rbit -x -f \$smaltoutput -o tmp1.\$smaltoutputsuffix ${meta.ID}_tmp.index ${name_1_fastq}
                cmdline="map -y ${params.nomapid}\$rbit -x -f \$smaltoutput -o tmp1.\$smaltoutputsuffix ${meta.ID}_tmp.index ${fastqdir}${meta.ID}.fastq"
            fi
        fi

        if [ "${newsmalt}" = "false" ]; then
            samtools view -b -S tmp1.sam -t ${smalt_indexes} > tmp1.bam
            rm tmp1.sam
        fi
    else
        cp ${bam} tmp1.bam
    fi
    """
}

process SMALT_INDEX {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/run-smalt:0.0.2'

    input:
    tuple val(meta), path(ref)

    output:
    tuple path(ref), path("${ref}.*")

    script:
    """
    if [ "${params.human}" == "True" ]; then
        smalt index -k 20 -s 13 ${meta.ID}_tmp.index ${ref}
    else 
        smalt index -k 13 -s 1 ${meta.ID}_tmp.index ${ref}
    fi
    samtools faidx ${ref}
    """
}
