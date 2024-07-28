include { RUN_BWA } from './../modules/BWA.nf'
include { RUN_SMALT } from './../modules/SMALT.nf'
include { RUN_SSAHA } from './../modules/RUN_SSAHA.nf'
include { MAKEPILEUP_FROM_SAM } from './../sub-workflows/MAKEPILEUP_FROM_SAM.nf'

workflow CALL_MAPPING {
    take:
    files
    tmpname
    ref
    ref_fai
    ref_index

    main:

    files = UNZIP_GZ(files, tmpname)

    files = UN_BAM (files)

    bwa_ch = files.combine(ref).combine(ref_fai)

    ref_index = ref_index.collect()

    if (params.program == "BWA") {
        files = RUN_BWA (bwa_ch, ref_index)
    } else if (params.program == "SMALT") {
        files = RUN_SMALT(tmpname, bwa_ch, ref_index)
    } else if (params.program == "SSAHA") {
        files = RUN_SSAHA(bwa_ch)
    }

    files = MAKEPILEUP_FROM_SAM(files)

}

process UNZIP_GZ {
    container 'zcat'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(files), path(file1), path(file2)
    val tmpname

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
    container 'bam_filter'
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