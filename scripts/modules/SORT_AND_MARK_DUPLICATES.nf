// include './../sub-workflows/makepileup_from_sam.nf'

process SORT_AND_MARK_DUPLICATES {
    container 'mpfs-gatk'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(pools), path(file1), path(file2), path(bam_file), val(cmdline)

    output:
    tuple val(pools), path(file1), path(file2), path("${runname}/tmp1.bam"), val(cmdline)
    path "${runname}/${name}_metrics.txt", emit: matrix_file_ch

    script:
    runname = pools.runname
    name = pools.name
    """
    mkdir -p ${runname}
    samtools sort ${bam_file} > ${runname}/tmpsort.bam
    #picard MarkDuplicates INPUT=${runname}/tmpsort.bam OUTPUT=${runname}/tmp1.bam METRICS_FILE=${runname}/${name}_metrics.txt
    gatk MarkDuplicates \
    -I ${runname}/tmpsort.bam \
    -O ${runname}/tmp1.bam \
    -M ${runname}/${name}_metrics.txt
    rm ${runname}/tmpsort.bam
    """
}
