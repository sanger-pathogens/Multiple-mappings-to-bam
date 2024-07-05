process INDEL_REALIGNMENT {
    input:
    tuple val(name), path (tmp1_bam)
    path ref

    output:
    tuple val(name), path ("${params.runname}/tmp1.bam")

    script:
    """
    mkdir -p "${params.runname}"
    samtools index ${tmp1_bam}
    cp ${ref} ${params.runname}/tmpref.fa
    samtools faidx ${params.runname}/tmpref.fa
    picard CreateSequenceDictionary R=${params.runname}/tmpref.fa O=${params.runname}/tmpref.dict
    gatk -I ${tmp1_bam} -R ${params.runname}/tmpref.fa -T RealignerTargetCreator -o ${params.runname}/tmp.intervals
    gatk -I ${tmp1_bam} -R ${params.runname}/tmpref.fa -T IndelRealigner --filter_bases_not_stored -targetIntervals ${params.runname}/tmp.intervals -o ${params.runname}/tmp.bam
    mv ${params.runname}/tmp.bam ${params.runname}/tmp1.bam
    rm *.bai ${params.runname}/tmpref.* ${params.runname}/tmp.intervals ${params.runname}/tmphead.*
    """
}