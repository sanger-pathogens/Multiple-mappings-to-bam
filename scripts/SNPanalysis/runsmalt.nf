params.smaltoutput="bam"
params.smaltoutputsuffix="bam"
params.rbit=""

process runSMALT {
    input:
        path bashfile
        val ref
        val fastqdir
        val name
        val pairedend
        val runname
        val domapping
        val newsmalt
        val tmpname
        val bam
    
    output:
        path bashfile

    script:
    """
    touch ${bashfile}
    if [ "${domapping}" = "true" ]; then
        echo "Running SMALT on ${name}..."

        if [ "${newsmalt}" = "true" ]; then
            ${params.smaltoutput}="bam"
            ${params.smaltoutputsuffix}="bam"
        else
            ${params.smaltoutput}="samsoft"
            ${params.smaltoutputsuffix}="sam"
        fi

        if [ "${pairedend}" = "true" ]; then
            if [ "${params.maprepeats}" = "true" ]; then
                echo "smalt map -y ${params.nomapid} -x -r 0 -i ${params.maxinsertsize} -j ${params.mininsertsize} -f ${params.smaltoutput} -o ${runname}/tmp1.${params.smaltoutputsuffix} ${tmpname}.index ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq" >> ${bashfile}
            else
                if [ "${newsmalt}" = "true" ]; then
                    ${params.rbit}=" -r -1"
                else
                    ${params.rbit}=""
                fi
                echo "smalt map -y ${params.nomapid}${params.rbit} -x -i ${params.maxinsertsize} -j ${params.mininsertsize} -f ${params.smaltoutput} -o ${runname}/tmp1.${params.smaltoutputsuffix} ${tmpname}.index ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq" >> ${bashfile}
            fi
        else
            if [ "${params.maprepeats}" = "true" ]; then
                echo "smalt map -y ${params.nomapid} -x -r 0 -f ${params.smaltoutput} -o ${runname}/tmp1.${params.smaltoutputsuffix} ${tmpname}.index ${fastqdir}${name}.fastq" >> ${bashfile}
            else
                if [ "${newsmalt}" = "true" ]; then
                    ${params.rbit}=" -r -1"
                else
                    ${params.rbit}=""
                fi
                echo "smalt map -y ${params.nomapid}${params.rbit} -x -f ${params.smaltoutput} -o ${runname}/tmp1.${params.smaltoutputsuffix} ${tmpname}.index ${fastqdir}${name}.fastq" >> ${bashfile}
            fi
        fi

        if [ "${newsmalt}" = "false" ]; then
            echo "samtools view -b -S ${runname}/tmp1.sam -t ${ref}.fai > ${runname}/tmp1.bam" >> ${bashfile}
            echo "rm ${runname}/tmp1.sam" >> ${bashfile}
        fi
    else
        echo "cp ${bam} ${runname}/tmp1.bam" >> ${bashfile}
    fi

    if [ "${fastqdir}" = "${tmpname}_unzipped/" ]; then
        echo "rm ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq" >> ${bashfile}
    fi
    """
}