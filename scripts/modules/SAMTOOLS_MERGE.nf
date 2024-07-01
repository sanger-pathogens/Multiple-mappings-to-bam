process SAMTOOLS_MERGE {
    input:
    val runname
    val name
    path name_bam
    path tmphead_sam

    output:
    path "${runname}/tmp1.bam", emit: tmp1_bam
    path "${runname}/tmphead.bam", emit; tmphead_bam

    script:
    """
    mkdir -p ${runname}
    samtools view -b -o ${runname}/tmphead.bam -H ${name_bam}
    samtools merge -c -p -f -r -h ${tmphead_sam} ${runname}/tmp.bam ${name_bam} ${runname}/tmphead.bam
    mv ${runname}/tmp.bam ${runname}/tmp1.bam
    rm ${runname}/${name}.bam
    """
}