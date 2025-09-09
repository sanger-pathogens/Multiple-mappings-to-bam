#!/usr/bin/env nextflow
// Copyright (C) 2024 Genome Research Ltd.

/*
========================================================================================
    HELP
========================================================================================
*/

def logo = NextflowTool.logo(workflow, params.monochrome_logs)

log.info logo

NextflowTool.commandLineParams(workflow.commandLine, log, params.monochrome_logs)


def printHelp() {
    NextflowTool.help_message("${workflow.ProjectDir}/schema.json", 
                               [],
    params.monochrome_logs, log)
}

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//FUNCTIONS
include { log_commandline           } from './modules/helper_functions.nf'

//MODULES
include { CONCAT_REFERENCE          } from './modules/concat_reference.nf'

//SUBWORKFLOWS
include { CALL_MAPPING              } from './sub-workflows/CALL_MAPPING.nf'
include { MAKE_PILEUP_FROM_SAM      } from './sub-workflows/MAKE_PILEUP_FROM_SAM.nf'
include { PSEUDOSEQUENCE_GENERATION } from './sub-workflows/PSEUDOSEQUENCE_GENERATION.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    if (params.help) {
        printHelp()
        exit 0
    }

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
    
    CALL_MAPPING(read_ch, params.ref)
    | MAKE_PILEUP_FROM_SAM
    | set { called_ch }

    if (params.pseudosequence == true) {
        PSEUDOSEQUENCE_GENERATION(called_ch, params.ref)
    }

}