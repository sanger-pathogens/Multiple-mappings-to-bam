process JOIN_DNA_INDELS {
    label "cpu_1"
    label "mem_16"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    container 'quay.io/ssd28/gsoc-experimental/join_dna_files_with_indels:0.0.2'

    input:
    path(mfa_list)
    path(ref)

    output:
    tuple path("${finalName}"), path(ref), emit: indel_joined_ch

    script:
    finalName="${ref.baseName}.aln"
    """
    join_dna_files_with_indels.py -r ${ref} -o ${finalName} -t ${mfa_list}
    """
}