include { PSEUDOSEQUENCE } from './../modules/PSEUDOSEQUENCE.nf'

workflow PSEUDOSEQUENCE_GENERATION {
    take:
    ref
    called_ch
    tmpname

    main:
    pseudosequence = PSEUDOSEQUENCE(called_ch)

    output_aln = JOIN_DNA_INDELS(pseudosequence, ref)

    SUMMARISE_SNPS(output_aln, ref)
}

process JOIN_DNA_INDELS {
    label "cpu_1"
    label "mem_16"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/join_dna_files_with_indels:0.0.2'

    input:
    tuple path(mfa), path(snp_filter_log), path(filtered_mapping_plot), path(indels)
    path ref

    output:
    path("${output}.aln")

    script:
    output = params.output
    if (params.indels == true) {
        """
        echo ${mfa} > mfa_list.txt
        join_dna_files_with_indels.py -r ${ref} -o ${output}.aln -t mfa_list.txt
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

    container 'quay.io/ssd28/gsoc-experimental/summarise_snps:0.0.2'

    input:
    path output_aln
    path ref

    script:
    summarystring = "summarise_snps.py -g -w -r "+ params.ref.split("/")[-1].split("\\.")[0] + " -o " + params.output + " -i " + params.output+".aln"
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