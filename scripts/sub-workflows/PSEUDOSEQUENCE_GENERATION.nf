include { PSEUDOSEQUENCE  } from './../modules/pseudosequence.nf'
include { JOIN_DNA_INDELS } from './../modules/join_dna_indels.nf'
include { SUMMARISE_SNPS  } from './../modules/summarise_snps.nf'

workflow PSEUDOSEQUENCE_GENERATION {
    
    take:
    ref
    called_ch

    main:
    pseudosequence = PSEUDOSEQUENCE(called_ch)

    output_aln = JOIN_DNA_INDELS(pseudosequence, ref)

    SUMMARISE_SNPS(output_aln, ref)
}