workflow PSEUDOSEQUENCE_GENERATION {
    take:
    files
    tmpname

    main:
    mfas_txt = MAKE_MFAS(tmpname)

    APPEND_MFA(files)




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
    publishDir "${params.outdir}", mode: 'copy'
    input:
    val file
    output:
    val 
    path("${tmpname}_mfas.txt")

    script:
}