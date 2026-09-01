process FILTER_BAM {
    tag "${meta.ID}"
    
    label "cpu_1"
    label "mem_1"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'

    publishDir "${params.outdir}/${meta.ID}_${params.program}/filtered_bams", mode: 'copy', overwrite: true, pattern: "${meta.ID}*.bam"

    input:
    tuple val(meta), path(tmp_bam)

    output:
    tuple val(meta), path("${meta.ID}.bam"), emit: bam_ch
    path "${meta.ID}_*.bam", optional: true

    script:
    if ( params.filter == 1 )
        """
        ln -s ${tmp_bam} ${meta.ID}.bam
        """
    else if( params.filter == 2 )
        """
        samtools view -F 4 -b -o ${meta.ID}.bam ${tmp_bam}
        """
    else if( params.filter == 3 )
        """
        samtools view -f 2 -b -o ${meta.ID}.bam ${tmp_bam}
        """
    else if( params.filter == 4 )
        """
        samtools view -F 4 -b -o ${meta.ID}.bam ${tmp_bam}
        samtools view -f 4 -b -o ${meta.ID}_unmapped.bam ${tmp_bam}
        """
    else if( params.filter == 5 ) 
        """
        samtools view -f 2 -b -o ${meta.ID}.bam ${tmp_bam}
        samtools view -F 2 -b -o ${meta.ID}_unpaired.bam ${tmp_bam}
        """
}