process INDEL_REALIGNMENT {
    input:
    path tmp1_bam
    path ref
    val runname

    output:
    path "${runname}/tmp1.bam"

    script:
    """
    mkdir -p "${runname}"
    samtools index ${tmp1_bam}
    cp ${ref} ${runname}/tmpref.fa
    samtools faidx ${runname}/tmpref.fa
    picard CreateSequenceDictionary R=${runname}/tmpref.fa O=${runname}/tmpref.dict
    gatk -I ${tmp1_bam} -R ${runname}/tmpref.fa -T RealignerTargetCreator -o ${runname}/tmp.intervals
    gatk -I ${tmp1_bam} -R ${runname}/tmpref.fa -T IndelRealigner --filter_bases_not_stored -targetIntervals ${runname}/tmp.intervals -o ${runname}/tmp.bam
    mv ${runname}/tmp.bam ${runname}/tmp1.bam
    rm *.bai ${runname}/tmpref.* ${runname}/tmp.intervals ${runname}/tmphead.*
    """
}