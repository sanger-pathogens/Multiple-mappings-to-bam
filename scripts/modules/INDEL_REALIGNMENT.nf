process INDEL_REALIGNMENT {
    label "cpu_1"
    label "mem_16"
    label "time_12"

    container 'quay.io/ssd28/gsoc-experimental/gatk-bcftools-samtools:0.0.2'
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pools), path(file1), path(file2), path (tmp1_bam), path(tmphead_sam), path(bam_bai), path (ref)

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
    
    gatk CreateSequenceDictionary \
    -R ${runname}/tmpref.fa \
    -O ${runname}/tmpref.dict

    # Generate GVCF file
    gatk HaplotypeCaller \
        -R ${runname}/tmpref.fa \
        -I ${tmp1_bam} \
        -O ${runname}/tmp.g.vcf \
        -bamout ${runname}/tmp.bam

    mv ${runname}/tmp.bam ${runname}/tmp1.bam
    """
}