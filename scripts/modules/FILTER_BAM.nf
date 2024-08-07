process FILTER_BAM {
    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/samtools:1.3'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pools), path(file1), path(file2), path (tmp_bam), path(tmphead_sam), path(bam_bai)

    output:
    tuple val(pools), path(file1), path(file2), path("${runname}/${name}.bam"), path(tmphead_sam), path(bam_bai)
    path "${runname}/${name}_*.bam", emit: filter_bam_ch, optional: true

    script:
    runname = pools.runname
    name = pools.name
    if ( params.filter == 1 )
        """
        mkdir -p "${runname}"
        mv ${tmp_bam} ${runname}/${name}.bam
        """
    else if( params.filter == 2 )
        """
        mkdir -p "${runname}"
        samtools view -F 4 -b -o ${runname}/${name}.bam ${tmp_bam}
        """
    else if( params.filter == 3 )
        """
        mkdir -p "${runname}"
        samtools view -f 2 -b -o ${runname}/${name}.bam ${tmp_bam}
        """
    else if( params.filter == 4 )
        """
        mkdir -p "${runname}"
        samtools view -F 4 -b -o ${runname}/${name}.bam ${tmp_bam}
        samtools view -f 4 -b -o ${runname}/${name}_unmapped.bam ${tmp_bam}
        """
    else if( params.filter == 5 ) 
        """
        mkdir -p "${runname}"
        samtools view -f 2 -b -o ${runname}/${name}.bam ${tmp_bam}
        samtools view -F 2 -b -o ${runname}/${name}_unpaired.bam ${tmp_bam}
        """
}