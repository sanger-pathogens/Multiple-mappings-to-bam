process INDEL_REALIGNMENT {
    tag "${meta.ID}"
    
    label "cpu_4"
    label "mem_10"
    label "time_12"

    container 'quay.io/ssd28/gsoc-experimental/gatk:3.7.0 '


    input:
    tuple val(meta), path(bam), path(bam_bai)
    tuple path(ref), path(fai), path(dict)

    output:
    tuple val(meta), path ("${meta.ID}_aligned.bam"), emit: bam_ch

    script:
    """
    java -jar /opt/GenomeAnalysisTK.jar -I ${bam}  -R ${ref} -T RealignerTargetCreator -o ${meta.ID}.intervals

    java -jar /opt/GenomeAnalysisTK.jar -I ${bam}  -R ${ref} -T IndelRealigner --filter_bases_not_stored -targetIntervals ${meta.ID}.intervals -o ${meta.ID}_aligned.bam
    """
}