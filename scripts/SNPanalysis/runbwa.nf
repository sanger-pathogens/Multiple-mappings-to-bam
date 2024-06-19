process RUN_BWA {
    input:
        path bashfile
        val fastqdir
        val name
        val pairedend
        val runname
        val domapping
        val is_zipped
        val bam
    
    output:
        path bashfile

    script:
    """
    # Create the bashfile if it does not exist
    touch ${bashfile}

    echo "Running BWA on ${name}..."

    if [ "${domapping}" = "true" ]; then
        if [ "${is_zipped}" = "true" ]; then
            if [ "${pairedend}" = "true" ]; then
                echo "bwa mem -v 1 -M -a -t 1 ${params.ref} ${fastqdir}${name}_1.fastq.gz ${fastqdir}${name}_2.fastq.gz > ${runname}/tmp.sam" >> ${bashfile}
            else
                echo "bwa mem -v 1 -M -a -t 1 ${params.ref} ${fastqdir}${name}.fastq.gz > ${runname}/tmp.sam" >> ${bashfile}
            fi
        else
            if [ "${pairedend}" = "true" ]; then
                echo "bwa mem -v 1 -M -a -t 1 ${params.ref} ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq > ${runname}/tmp.sam" >> ${bashfile}
            else
                echo "bwa mem -v 1 -M -a -t 1 ${params.ref} ${fastqdir}${name}.fastq > ${runname}/tmp.sam" >> ${bashfile}
            fi
        fi

        echo "samtools view -b -S ${runname}/tmp.sam -t ${params.ref}.fai > ${runname}/tmp1.bam" >> ${bashfile}
        echo "rm -f ${runname}/tmp.sam" >> ${bashfile}
    else
        echo "cp ${bam} ${runname}/tmp1.bam" >> ${bashfile}
    fi
    """
}