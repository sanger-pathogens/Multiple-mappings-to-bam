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

    if (!params.read_dir) {
        exit 1, 'Error: Please provide a read directory using --read_dir'
    }

    Channel.fromFilePairs("${params.read_dir}/*_{1,2}.fastq.gz")
    | map { id, reads ->
        meta = [:]
        meta.ID = id
        [meta, reads[0], reads[1]]
    }
    | set { read_ch }

    Channel.fromPath(params.ref)
    | set { reference_ch }

    if (params.cat_reference) {
        CONCAT_REFERENCE(reference_ch)
        | set{ ref }
    } else {

        reference_ch
        | set { ref }
    }
    
    CALL_MAPPING(read_ch, ref)
    | MAKE_PILEUP_FROM_SAM
    | set { called_ch }

    if (params.pseudosequence == true) {
        PSEUDOSEQUENCE_GENERATION(called_ch, ref)
    }

}