// include './../sub-workflows/makepileup_from_sam.nf'

process SORT_AND_MARK_DUPLICATES {
    // container 'quay.io/ssd28/gsoc-experimental/sort-and-mark-duplicates:0.0.1'

    input:
        tuple val(name), path(bam_file)

    output:
        tuple val(name), path ("${params.runname}/tmp1.bam"), emit: tmp1_bam
        path "${params.runname}/${name}_metrics.txt", emit: matrix_file_ch

    script:
    """
    mkdir -p ${params.runname}
    samtools sort ${bam_file} > ${params.runname}/tmpsort.bam
    picard MarkDuplicates INPUT=${params.runname}/tmpsort.bam OUTPUT=${params.runname}/tmp1.bam METRICS_FILE=${params.runname}/${name}_metrics.txt
    rm ${params.runname}/tmpsort.bam
    """
}
