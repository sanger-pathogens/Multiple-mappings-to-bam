process INDEL_REALIGNMENT {
    container 'mpfs-gatk'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pools), path(file1), path(file2), path (tmp1_bam), path(tmphead_sam), path(bam_bai)
    path ref

    output:
    tuple val(pools), path(file1), path(file2), path ("${runname}/tmp1.bam"), path(tmphead_sam), path(bam_bai)

    script:
    runname = pools.runname
    name = pools.name
    """
    mkdir -p "${runname}"
    samtools index ${tmp1_bam}
    cp ${ref} ${runname}/tmpref.fa
    samtools faidx ${runname}/tmpref.fa

    #picard CreateSequenceDictionary R=${runname}/tmpref.fa O=${runname}/tmpref.dict
    gatk CreateSequenceDictionary \
    -R ${runname}/tmpref.fa \
    -O ${runname}/tmpref.dict

    #gatk -I ${tmp1_bam} -R ${runname}/tmpref.fa -T RealignerTargetCreator -o ${runname}/tmp.intervals
    #gatk -I ${tmp1_bam} -R ${runname}/tmpref.fa -T IndelRealigner --filter_bases_not_stored -targetIntervals ${runname}/tmp.intervals -o ${runname}/tmp.bam

    # Generate GVCF file
    gatk HaplotypeCaller \
        -R ${runname}/tmpref.fa \
        -I ${tmp1_bam} \
        -O ${runname}/tmp.g.vcf \
        -bamout ${runname}/tmp.bam


    mv ${runname}/tmp.bam ${runname}/tmp1.bam
    rm *.bai ${runname}/tmpref.* ${runname}/tmp.intervals ${runname}/tmphead.*
    """
}