process RUN_BWA {
    container 'quay.io/ssd28/gsoc-experimental/run-bwa:0.0.1'
    tag "${name}"
    input:
        val fastqdir
        val name
        val pairedend
        val runname
        val domapping
        val is_zipped
        val bam

    script:
    """

    #if [ "${domapping}" = "true" ]; then
        if [ "${is_zipped}" = "true" ]; then
            if [ "${pairedend}" = "true" ]; then
                bwa mem -v 1 -M -a -t 1 ${params.ref} ${fastqdir}${name}_1.fastq.gz ${fastqdir}${name}_2.fastq.gz > ${runname}/tmp.sam
            else
                bwa mem -v 1 -M -a -t 1 ${params.ref} ${fastqdir}${name}.fastq.gz > ${runname}/tmp.sam
            fi
        else
            if [ "${pairedend}" = "true" ]; then
                bwa mem -v 1 -M -a -t 1 ${params.ref} ${fastqdir}${name}_1.fastq ${fastqdir}${name}_2.fastq > ${runname}/tmp.sam
            else
                bwa mem -v 1 -M -a -t 1 ${params.ref} ${fastqdir}${name}.fastq > ${runname}/tmp.sam
            fi
        fi

        samtools view -b -S ${runname}/tmp.sam -t ${params.ref}.fai > ${runname}/tmp1.bam
        rm -f ${runname}/tmp.sam
    #else
    #    cp ${bam} ${runname}/tmp1.bam
    #fi
    """
}
