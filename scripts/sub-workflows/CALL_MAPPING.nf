include { RUN_BWA } from './../modules/RUN_BWA.nf'

workflow CALL_MAPPING {
    take:
    files
    ziplist
    tmpname

    main:

    gz_files = ziplist
        .map { it.values() }
        .flatten()

    bam_ch = files.map { row ->
        def bam_path = ""
        if (row[2] == "") {
            bam_path = baseDir+"/"+params.nfile
        } else {
            bam_path = row[2]
        }
        [bam_path, row[0].fastqdir, row[0].name]
    }
    
    UNZIP_GZ(gz_files, fastqdir_zip)

    UN_BAM (bam_ch)
    
}

process UNZIP_GZ {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gz_file
    val fastqdir

    output:
    path "${outputFilename}"

    when:
    params.keep == false

    script:
    outputFilename = fastqdir + gz_file.getName().split('\\.')[0..-2].join('.')
    """
    mkdir -p ${fastqdir}
    zcat ${gz_file} > ${outputFilename}
    """
}

process UN_BAM {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(bamlist), val(fastqdir), val(name)

    output:
    path "${outputFilename}"

    when:
    params.keep == false && bamlist != ""

    script:
    outputFilename = "${fastqdir}${name}"
    
    if (params.pairedend) {
        """
        mkdir -p ${fastqdir}
        bam_filter.py -t all -b ${bamlist} -o ${outputFilename}
        """
    } else {
        """
        mkdir -p ${fastqdir}
        bam_filter.py -t all -f fasta -b ${bamlist} -o ${outputFilename}
        """
    }
}