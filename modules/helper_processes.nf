process UNZIP_GZ {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_100M"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/zcat:0.0.2'

    input:
    tuple val(meta), path(file1), path(file2)

    output:
    tuple val(meta), path(outputFileName), path(outputFileName2)

    script:
    outputFileName = file1.name.endsWith('.gz') ? file1.baseName : file1
    outputFileName2 = file2.name.endsWith('.gz') ? file2.baseName : file2

    """
    gunzip -f ${file1} || true
    gunzip -f ${file2} || true
    """
}

process UN_BAM {
    label "cpu_1"
    label "mem_100M"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/bam_filter:0.0.2'
    
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(outputFileName1), path(outputFileName2)

    script:
    outputFileName1 = "${meta.ID}_1.fastq"
    outputFileName2 = "${meta.ID}_2.fastq"

    if (params.domapping == true && bam.name.endsWith('.bam')) {
        """
        bam_filter.py -t all -b ${bam} -o ${meta.ID}
        """
    } else {
        """
        """
    }
}