process UNZIP {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(outdir), val(files)

    script:

    options = files[0]
    ziplist = files[1]
    println outdir
    println files

    unZipCmd = ["echo 'unzipping'"]

    is_zip = ziplist.containsKey(options.name)

    if (is_zip) {
        ziplist[options.name].each { filename ->
            outputFilename = options.fastqdir + filename.split('/')[-1].split('\\.')[0..-2].join('.')
            command = "zcat ${filename} > ${outputFilename}"
            unZipCmd << command
        }
    }

    bam_exist = file("${outdir}/${options.runname}/${options.name}.bam").exists()
    bcf_exist = file("${outdir}/${options.runname}/${options.name}.bcf").exists()
    mfa_exist = file("${outdir}/${options.runname}/${options.name}.mfa").exists()
    
    if (params.keep && params.pseudosequence && bam_exist && bcf_exist && mfa_exist) {
        """
        echo "${options.runname} contains a .bam, .bcf and .mfa file already. As you selected the keep option, it will not be mapped again." 
        """
    } else if (params.keep && !params.pseudosequence && bam_exist && bcf_exist) {
        """
        echo "${options.runname} contains a .bam, and .bcf file already. As you selected the keep option, it will not be mapped again."
        """
    } else {
        """
        if [ "${is_zip}" == "true" ]; then
            ${unZipCmd.join('\n')}
        else 
            echo "No file found"
        fi
        """
    }
}