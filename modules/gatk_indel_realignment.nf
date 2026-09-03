// GATK requires a readgroup (@RG) tag for indel realignment, so we add one to the BAM header if it doesn't already have one
process ADD_READGROUP {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_100M"
    label "time_30m"

    container 'quay.io/ssd28/gsoc-experimental/void:0.0.1'

    stageInMode 'copy'

    input:
    tuple val(meta), path(header)
    val(program)

    output:
    tuple val(meta), path(header), emit: header_ch

    script:
    """
    if ! grep -q "^@RG" "${header}"; then
        now=\$(date +'%Y-%m-%dT%H:%M:%S')
        echo '@RG\tID:${meta.ID}\tCN:Sanger\tDT:'\$now'\tPG:${program}\tPL:ILLUMINA\tSM:${meta.ID}' >> ${header}
    fi
    """
}

process INDEL_REALIGNMENT {
    tag "${meta.ID}"
    
    label "cpu_4"
    label "mem_10"
    label "time_12"

    container 'quay.io/ssd28/gsoc-experimental/gatk:3.7.0'

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