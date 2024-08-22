process FORMAT_SMALT_HEADER {
    label "cpu_1"
    label "mem_16"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/void:0.0.1'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pools), path(file1), path(file2), path (name_bam), val(cmdline), path(tmphead_sam), path(bam_bai)

    output:
    tuple val(pools), path(file1), path(file2), path (name_bam), path(tmphead_sam), path(bam_bai)

    script:
    name = pools.name
    newsmalt = "true"
    """
    now=\$(date +'%Y-%m-%dT%H:%M:%S')
    echo "@RG\tID:${name}\tCN:Sanger\tDT:"\$now"\tPG:SMALT\tPL:ILLUMINA\tSM:${name}" >> ${tmphead_sam}
    if [ ${params.domapping} ] && [ ${newsmalt} = "false" ]
    then
        smaltversion=\$( smalt version | grep Version | awk '{print \$2}' )
        echo "@PG\tID:SMALT\tPN:SMALT\tCL:${cmdline}\tVN:\$smaltversion" >> ${tmphead_sam}
    fi
    """
}


process RUN_SMALT {
    label "cpu_1"
    label "mem_16"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/run-smalt:0.0.2'

    tag "${name}"
    
    input:
    tuple val(pools), path(name_1_fastq), path(name_2_fastq), path(ref), path(ref_fai)
    tuple path(f1), path(f2)

    output:
    tuple val(pools), path(name_1_fastq), path(name_2_fastq), path("${runname}/tmp1.bam"), env(cmdline)


    script:
    newsmalt = true
    runname = pools.runname
    domapping = pools.domapping
    pairedend = pools.pairedend
    fastqdir = pools.fastqdir
    name = pools.name
    bam = pools.bam
    """
    mkdir -p "${runname}"

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
                smalt map -y ${params.nomapid} -x -r 0 -i ${params.maxinsertsize} -j ${params.mininsertsize} -f \$smaltoutput -o ${runname}/tmp1.\$smaltoutputsuffix ${name}_tmp.index ${name_1_fastq} ${name_2_fastq}
                cmdline="map -y ${params.nomapid} -x -r 0 -i ${params.maxinsertsize} -j ${params.mininsertsize} -f \$smaltoutput -o ${runname}/tmp1.\$smaltoutputsuffix ${name}_tmp.index ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq"
            else
                if [ "${newsmalt}" = "true" ]; then
                    rbit=" -r -1"
                else
                    rbit=""
                fi
                smalt map -y ${params.nomapid}\$rbit -x -i ${params.maxinsertsize} -j ${params.mininsertsize} -f \$smaltoutput -o ${runname}/tmp1.\$smaltoutputsuffix ${name}_tmp.index ${name_1_fastq} ${name_2_fastq}
                cmdline="map -y ${params.nomapid}\$rbit -x -i ${params.maxinsertsize} -j ${params.mininsertsize} -f \$smaltoutput -o ${runname}/tmp1.\$smaltoutputsuffix ${name}_tmp.index ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq"
            fi
        else
            if [ "${params.maprepeats}" = "true" ]; then
                smalt map -y ${params.nomapid} -x -r 0 -f \$smaltoutput -o ${runname}/tmp1.\$smaltoutputsuffix ${name}_tmp.index ${name_1_fastq}
                cmdline="map -y ${params.nomapid} -x -r 0 -f \$smaltoutput -o ${runname}/tmp1.\$smaltoutputsuffix ${name}_tmp.index ${fastqdir}${name}.fastq"
            else
                if [ "${newsmalt}" = "true" ]; then
                    \$rbit=" -r -1"
                else
                    \$rbit=""
                fi
                smalt map -y ${params.nomapid}\$rbit -x -f \$smaltoutput -o ${runname}/tmp1.\$smaltoutputsuffix ${name}_tmp.index ${name_1_fastq}
                cmdline="map -y ${params.nomapid}\$rbit -x -f \$smaltoutput -o ${runname}/tmp1.\$smaltoutputsuffix ${name}_tmp.index ${fastqdir}${name}.fastq"
            fi
        fi

        if [ "${newsmalt}" = "false" ]; then
            samtools view -b -S ${runname}/tmp1.sam -t ${ref_fai} > ${runname}/tmp1.bam
            rm ${runname}/tmp1.sam
        fi
    else
        cp ${bam} ${runname}/tmp1.bam
    fi
    """
}

process SMALT_INDEX {
    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/run-smalt:0.0.2'

    input:
    path ref

    output:
    tuple path("*.smi"), path("*.sma")
    path("*.fai"), emit: fai

    script:
    """
    if [ "${params.human}" == "True" ]; then
        smalt index -k 20 -s 13 ${name}_tmp.index ${ref}
    else 
        smalt index -k 13 -s 1 ${name}_tmp.index ${ref}
    fi
    samtools faidx ${ref}
    """
}