include { PSEUDOSEQUENCE  } from './../modules/pseudosequence.nf'
include { JOIN_DNA_INDELS } from './../modules/join_dna_indels.nf'
include { SUMMARISE_SNPS  } from './../modules/summarise_snps.nf'

workflow PSEUDOSEQUENCE_GENERATION {
    
    take:
    called_ch
    ref

    main:
    PSEUDOSEQUENCE(called_ch)

    JOIN_DNA_INDELS(PSEUDOSEQUENCE.out.pseudosequence, ref)
    | SUMMARISE_SNPS

}