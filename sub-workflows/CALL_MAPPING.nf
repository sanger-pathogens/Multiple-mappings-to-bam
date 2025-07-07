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

    // Filepath ref is a string, convert to file object
    ref = file(ref)

    switch (params.program.toUpperCase()) {
        case "BWA":
            BWA_INDEX(ref)
            | set { ref_plus_index }

            RUN_BWA(unzipped_reads, ref_plus_index)
            | set { mapped_ch }
            
            break

        case "SMALT":
            SMALT_INDEX(ref)
            | set { ref_plus_index }

            RUN_SMALT(unzipped_reads, ref_plus_index)
            | set { mapped_ch }

            break

        case "SSAHA":
            reads_and_ref_ch = unzipped_reads.combine(ref).combine(index_ch)
            mapped_ch = RUN_SSAHA(reads_and_ref_ch)

            break

        default:
            log.error("Unsupported program: ${params.program}")
    }

    emit:
    mapped_ch
    ref_plus_index
}
