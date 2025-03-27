//FUNCTIONS
include { log_commandline           } from './modules/helper_functions.nf'

//MODULES
include { CONCAT_REFERENCE          } from './modules/concat_reference.nf'

//SUBWORKFLOWS
include { CALL_MAPPING              } from './sub-workflows/CALL_MAPPING.nf'
include { MAKE_PILEUP_FROM_SAM      } from './sub-workflows/MAKE_PILEUP_FROM_SAM.nf'
include { PSEUDOSEQUENCE_GENERATION } from './sub-workflows/PSEUDOSEQUENCE_GENERATION.nf'

workflow {
    log_commandline()

    Channel.fromFilePairs("/lustre/scratch126/pam/teams/team230/sd28/Multiple-mappings-to-bam/test_cases/reads/*_{1,2}.fastq.gz")
    | map { id, reads ->
        meta = [:]
        meta.ID = id
        [meta, reads[0], reads[1]]
    }
    | set { read_ch }

    CONCAT_REFERENCE(params.ref)
    | set{ ref }
    
    CALL_MAPPING(read_ch, params.ref)
    | MAKE_PILEUP_FROM_SAM
    | set { called_ch }

    if (params.pseudosequence == true) {
        PSEUDOSEQUENCE_GENERATION(params.ref, called_ch)
    }

}