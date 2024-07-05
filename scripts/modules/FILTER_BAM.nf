process FILTER_BAM {
    input:
    tuple val(name), path (tmp_bam)

    output:
    tuple val(name), path ("${params.runname}/${name}.bam"), emit: name_bam
    path "${params.runname}/${name}_*.bam", emit: filter_bam_ch, optional: true

    script:

    if ( params.filter == '1' )
        """
        mkdir -p "${params.runname}"
        mv ${tmp_bam} ${params.runname}/${name}.bam
        """
    else if( params.filter == '2' )
        """
        mkdir -p "${params.runname}"
        samtools view -F 4 -b -o ${params.runname}/${name}.bam ${tmp_bam}
        """
    else if( params.filter == '3' )
        """
        mkdir -p "${params.runname}"
        samtools view -f 2 -b -o ${params.runname}/${name}.bam ${tmp_bam}
        """
    else if( params.filter == '4' )
        """
        mkdir -p "${params.runname}"
        samtools view -F 4 -b -o ${params.runname}/${name}.bam ${tmp_bam}
        samtools view -f 4 -b -o ${params.runname}/${name}_unmapped.bam ${tmp_bam}
        """
    else if( params.filter == '5' ) 
        """
        mkdir -p "${params.runname}"
        samtools view -f 2 -b -o ${params.runname}/${name}.bam ${tmp_bam}
        samtools view -F 2 -b -o ${params.runname}/${name}_unpaired.bam ${tmp_bam}
        """
}