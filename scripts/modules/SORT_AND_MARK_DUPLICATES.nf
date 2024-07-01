// include './../sub-workflows/makepileup_from_sam.nf'

process SORT_AND_MARK_DUPLICATES {
    // container 'quay.io/ssd28/gsoc-experimental/sort-and-mark-duplicates:0.0.1'

    input:
        val runname
        val name
        path bam_file //${runname}/tmp1.bam

    output:
        path "${runname}/tmp1.bam", emit: tmp1_bam
        path "${runname}/${name}_metrics.txt", emit: matrix_file_ch

    script:
    """
    mkdir -p ${runname}
    samtools sort ${bam_file} > ${runname}/tmpsort.bam
    picard MarkDuplicates INPUT=${runname}/tmpsort.bam OUTPUT=${runname}/tmp1.bam METRICS_FILE=${runname}/${name}_metrics.txt
    rm ${runname}/tmpsort.bam
    """
}
