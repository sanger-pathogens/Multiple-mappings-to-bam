include { UNZIP_GZ
          UN_BAM                 } from './../modules/helper_processes.nf'
include { BWA_INDEX; RUN_BWA     } from './../modules/bwa.nf'
include { SMALT_INDEX; RUN_SMALT } from './../modules/smalt.nf'
include { RUN_SSAHA              } from './../modules/ssaha.nf'

workflow CALL_MAPPING {
    take:
    read_ch
    ref

    main:
    UNZIP_GZ(read_ch)
    | set { unzipped_reads }


    switch (params.program.toUpperCase()) {
        case "BWA":
            BWA_INDEX(ref)
            | set { ref_plus_index }

            RUN_BWA(unzipped_reads, ref_plus_index)
            | set { mapped_ch }
            
            break

        case "SMALT":
            (index_ch, fai) = SMALT_INDEX(ref)

            mapped_ch = RUN_SMALT(reads_and_ref_ch)
            break

        case "SSAHA":
            mapped_ch = RUN_SSAHA(reads_and_ref_ch)

            break

        default:
            log.error("Unsupported program: ${params.program}")
    }

    emit:
    mapped_ch
    ref_plus_index
}