include { RUN_BWA } from './../modules/BWA.nf'
include { RUN_SMALT } from './../modules/SMALT.nf'
include { RUN_SSAHA } from './../modules/RUN_SSAHA.nf'

workflow CALL_MAPPING {
    take:
    read_ch
    ref
    ref_fai
    ref_index

    main:

    unzipped_ch = UNZIP_GZ(read_ch)

    unbammed_ch = UN_BAM (unzipped_ch)

    reads_and_ref_ch = unbammed_ch.combine(ref).combine(ref_fai)

    ref_index = ref_index.collect()

    if (params.program == "BWA") {
        mapped_ch = RUN_BWA (reads_and_ref_ch, ref_index)
    } else if (params.program == "SMALT") {
        mapped_ch = RUN_SMALT(reads_and_ref_ch, ref_index)
    } else if (params.program == "SSAHA") {
        mapped_ch = RUN_SSAHA(reads_and_ref_ch)
    }

    emit:
    mapped_ch

}

process UNZIP_GZ {
    label "cpu_1"
    label "mem_1"
    label "time_1"
    
    container 'quay.io/ssd28/gsoc-experimental/zcat:0.0.2'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(files), path(file1), path(file2)

    output:
    tuple val(files), path(outputFileName), path(outputFileName2)

    script:
    fastqdir = files.fastqdir
    outputFileName=file1
    outputFileName2=file2
    
    paramsNfile = params.nfile
    if (params.keep == false && params.program != "BWA" && file1.name.endsWith('.fastq.gz')) {
        outputFileName = fastqdir+file1.name.split('\\.')[0..-2].join('.')
        outputFileName2 = fastqdir+file2.name.split('\\.')[0..-2].join('.')
        """
        mkdir -p ${fastqdir}
        zcat ${file1} > ${outputFileName}

        if [ "${file2.name}" != "${paramsNfile}" ]; then
            zcat ${file2} > ${outputFileName2}
        else
            cp ${file2} ${outputFileName2}
        fi
        """
    } else {
        """
        """
    }
}

process UN_BAM {
    label "cpu_1"
    label "mem_16"
    label "time_1"

    container 'quay.io/ssd28/gsoc-experimental/bam_filter:0.0.2'
    
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(pools), path(file1), path(file2)

    output:
    tuple val(pools), path(outputFileName), path(outputFileName2)

    script:
    outputFileName = file1
    outputFileName2 = file2
    fastqdir = pools.fastqdir
    name = pools.name
    
    if (!(params.keep == false && params.program != "BWA" && file1.name.endsWith('.fastq.gz')) && 
        !(params.program == "BWA" && file1.name.endsWith('.gz')) && 
        file1.name.endsWith('.bam') && 
        params.domapping == true) {
        outputFileName = "${fastqdir}${name}"
        if (params.pairedend) {
        """
        mkdir -p ${fastqdir}
        bam_filter.py -t all -b ${file1} -o ${outputFileName}
        """
        } else {
            """
            mkdir -p ${fastqdir}
            bam_filter.py -t all -f fasta -b ${file1} -o ${outputFileName}
            """
        }
    } else {
        """
        """
    }
    
    
}