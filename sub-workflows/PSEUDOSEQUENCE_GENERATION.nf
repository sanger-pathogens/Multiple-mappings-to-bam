include { PSEUDOSEQUENCE  } from './../modules/pseudosequence.nf'
include { JOIN_DNA_INDELS } from './../modules/join_dna_indels.nf'
include { SUMMARISE_SNPS  } from './../modules/summarise_snps.nf'

workflow PSEUDOSEQUENCE_GENERATION {
    
    take:
    called_ch
    ref

    main:
    PSEUDOSEQUENCE(called_ch)

    PSEUDOSEQUENCE.out.pseudosequence
    | map { meta, file -> file.toString() } //store the file name not the contents
    | collectFile(name: 'mfa_list.txt', newLine: true)
    | set { mfa_list }

    JOIN_DNA_INDELS(mfa_list, ref)
    | SUMMARISE_SNPS

}