workflow PSEUDOSEQUENCE_GENERATION {
    take:
    files
    tmpname

    main:
    ref = Channel.fromPath(params.ref)
    mfas_txt = MAKE_MFAS(tmpname)
    output_ch = files.map { row -> 
        row[0].runname+'/'+row[0].name+'_mfas.txt'
    }
    ls_ch=output_ch.toList()
    mfas_txt = APPEND_MFA(ls_ch, mfas_txt)

    ls_space = output_ch.toList()

    output_aln = JOIN_DNA_INDELS(mfas_txt, ref, ls_space)

    SUMMARISE_SNPS(output_aln, ref)
}

process MAKE_MFAS {
    input:
    val tmpname

    output:
    path("${tmpname}_mfas.txt")

    script:
    """
    touch ${tmpname}_mfas.txt
    """
}
process APPEND_MFA {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    val files
    path txt

    output:
    path(txt)

    script:
    str = files.join('\n')
    """
    echo "${str}" >> ${txt}
    """
}

process JOIN_DNA_INDELS {
    publishDir "${params.outdir}", mode: 'copy'
    container 'join_dna_files_with_indels'
    input:
    path (mfas_txt)
    path ref
    val ls_space

    output:
    path("${output}.aln")

    script:
    output = params.output
    str = ls_space.join(' ')
    if (params.indels == true) {
        """
        python2 /opt/join_dna_files_with_indels/join_dna_files_with_indels.py -r ${ref} -o ${output}.aln -t ${mfas_txt}
        """
    } else if (params.incref == false) {
        """
        cat ${str} > ${output}.aln
        """
    } else {
        """
        cat ${ref} ${str} > ${output}.aln
        """
    }
}

process SUMMARISE_SNPS {
    container 'summarise_snps'
    input:
    path output_aln
    path ref

    script:
    summarystring = "python2 /opt/summarise_snps/summarise_snps.py -g -w -r "+ params.ref.split("/")[-1].split("\\.")[0] + " -o " + params.output + " -i " + params.output+".aln"
    bootstrap = params.bootstrap

    if (params.embl != "") {
        summarystring = summarystring + " -e "+params.embl
    }
    if (params.alnfile == true) {
        summarystring = summarystring + " -a"
    }
    if (params.tabfile == true) {
        summarystring = summarystring + " -t"
    }
    if (params.raxml == true) {
        summarystring = summarystring + " -p -l -b ${bootstrap}"
    }

    """
    ${summarystring}
    """
}