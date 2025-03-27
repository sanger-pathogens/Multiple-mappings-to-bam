process RUN_SSAHA {
    tag "${meta.ID}"

    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    container 'quay.io/sangerpathogens/ssaha2:v2.5.5_cv3'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(name_1_fastq), path(name_2_fastq), path(ref), path(ref_fai)

    output:
    tuple val(meta), path(name_1_fastq), path(name_2_fastq), path("tmp1.bam")

    script:
    runname = meta.runname
    pairedend = meta.pairedend
    cmdline=""
    """
    if [ "${pairedend}" = "false" ]; then
        ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -output sam_soft -outfile tmp1.sam ${ref} ${name_fastq}
    else
        ssaha2 -score ${params.ssahaquality} -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -outfile tmp1.sam -pair ${params.mininsertsize},${params.maxinsertsize} -output sam_soft ${ref} ${name_1_fastq} ${name_2_fastq}
    fi

    samtools view -b -S tmp1.sam -t ${ref_fai} > tmp1.bam

    if [ "${pairedend}" = "true" ] && [ "${params.circular}" = "true" ]; then
        fix_circular_bams.py -b tmp1.bam -o tmp
        rm tmp1.bam
    else
        mv tmp1.bam tmp.bam
    fi

    samtools view -H tmp.bam > tmp2.sam
    cat tmp2.sam tmp1.sam > tmp.sam
    samtools view -b -S tmp.sam -t ${ref_fai} > tmp1.bam
    rm -f tmp.sam
    rm tmp2.sam tmp1.sam
    """
}
