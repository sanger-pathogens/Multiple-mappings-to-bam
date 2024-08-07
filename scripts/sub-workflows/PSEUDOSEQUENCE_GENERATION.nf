workflow PSEUDOSEQUENCE_GENERATION {
    take:
    files
    tmpname

    main:
    ref = Channel.fromPath(params.ref)

    mfas_txt = MAKE_TXT(tmpname)

    output_ch = files.map { row -> 
        row[0].name+'_mfas.mfa'
    }

    (mfa_files, txt_files) = MAKE_MFA(files)

    mfa_files = mfa_files.collect()
    txt_files = txt_files.collect()
    ls_ch=output_ch.toList()
    
    mfas_txt = APPEND_MFA(ls_ch, mfas_txt)

    ls_space = output_ch.toList()

    output_aln = JOIN_DNA_INDELS(mfas_txt, ref, mfa_files, ls_space, txt_files)

    SUMMARISE_SNPS(output_aln, ref)
}

process MAKE_MFA {
    label "cpu_1"
    label "mem_16"
    label "time_1"

    input:
    val files

    output:
    path "${runname}/${name}_mfas.mfa"
    path "${runname}/${name}_mfas_indels.txt"

    script:
    runname = files[0].runname
    name = files[0].name
    """
    mkdir -p ${runname}
    touch ${runname}/${name}_mfas.mfa
    touch ${runname}/${name}_mfas_indels.txt
    """
}
process MAKE_TXT {
    label "cpu_1"
    label "mem_16"
    label "time_1"

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
    label "cpu_1"
    label "mem_16"
    label "time_1"

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
    label "cpu_1"
    label "mem_16"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/join_dna_files_with_indels:0.0.1'

    input:
    path (mfas_txt)
    path ref
    path mfa_files
    val ls_space
    path txt_files

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
    label "cpu_1"
    label "mem_16"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/summarise_snps:0.0.1'

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