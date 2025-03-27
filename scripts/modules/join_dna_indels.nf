process JOIN_DNA_INDELS {
    label "cpu_1"
    label "mem_16"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/join_dna_files_with_indels:0.0.2'

    input:
    tuple val(meta), path(mfa)
    path(ref)

    output:
    path("${meta.ID}.aln"), path(ref), emit: indel_joined_ch

    script:
    if (params.indels == true) {
        """
        echo ${mfa} > mfa_list.txt
        join_dna_files_with_indels.py -r ${ref} -o ${meta.ID}.aln -t mfa_list.txt
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