include { UNZIP_GZ                                      } from '../modules/helper_processes.nf'
include { BWA_INDEX; RUN_BWA                            } from '../modules/bwa.nf'
include { SMALT_INDEX; 
          RUN_SMALT                                     } from '../modules/smalt.nf'
include { SAMTOOLS_SAM_TO_BAM as SAMTOOLS_SMALT_TO_BAM;
          SAMTOOLS_SAM_TO_BAM as SAMTOOLS_SSAHA_TO_BAM;
          INDEX_REF                                     } from '../modules/samtools.nf'
include { RUN_SSAHA;
          FIX_CIRCULAR_BAMS                             } from '../modules/ssaha.nf'

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
            SMALT_INDEX(ref)
            | set { ref_plus_index }

            RUN_SMALT(unzipped_reads, ref_plus_index)
            | SAMTOOLS_SMALT_TO_BAM
            | set { mapped_ch }

            break

        case "SSAHA":
            INDEX_REF(ref)
            | set { ref_plus_index }

            RUN_SSAHA(unzipped_reads, ref_plus_index)
            | SAMTOOLS_SSAHA_TO_BAM
            | set { raw_mapped_ch }

            if (params.circular) {

                FIX_CIRCULAR_BAMS(raw_mapped_ch)
                | set { mapped_ch }

            } else {

                raw_mapped_ch
                | set { mapped_ch }

            }

            break

        default:
            log.error("Unsupported program: ${params.program}")
    }

    emit:
    mapped_ch
    ref
}