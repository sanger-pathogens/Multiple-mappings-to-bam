include { MARK_DUPLICATES; SEQUENCE_DICT                } from '../modules/picard.nf'
include { INDEX_REF
          SAMTOOLS_SORT;
          SAMTOOLS_SORT as SAMTOOLS_SORT_RAW_MAPPING;
          UPDATE_HEADER_AND_RG_TAG;
          SAMTOOLS_INDEX as SAMTOOLS_RAW_INDEX;
          SAMTOOLS_INDEX as SAMTOOLS_FILTERED_INDEX;
          SAMTOOLS_PILEUP                               } from '../modules/samtools.nf'
include { INDEL_REALIGNMENT;                            } from '../modules/gatk_indel_realignment.nf'
include { FILTER_BAM                                    } from '../modules/filter_bam.nf'
include { BCFTOOLS_CALL                                 } from '../modules/bcftools.nf'

workflow MAKE_PILEUP_FROM_SAM {

    take:
    mapped_sam_ch
    ref

    main:

    SAMTOOLS_SORT_RAW_MAPPING(mapped_sam_ch)

    if (params.markdup) {

        MARK_DUPLICATES(SAMTOOLS_SORT_RAW_MAPPING.out.bam_ch)

        MARK_DUPLICATES.out.deduped_ch
        | set { deduped_ch }

    } else {

        SAMTOOLS_SORT_RAW_MAPPING.out.bam_ch
        | set{ deduped_ch }

    }

    // Map to generate program label PG for readgroup (RG)
    Map<String, String> program_label = [
        "BWA": "BWA MEM",
        "SSAHA": "SSAHA",
        "SMALT": "SMALT"
    ]

    UPDATE_HEADER_AND_RG_TAG(
        deduped_ch,
        Channel.value(program_label[params.program])
    )
    | SAMTOOLS_RAW_INDEX
    | set { sam_ref_ch }

    if (params.GATK) {

        INDEX_REF(ref)
        | SEQUENCE_DICT

        INDEL_REALIGNMENT(sam_ref_ch, SEQUENCE_DICT.out.ref_ch)
        | set { indel_realigned_ch }

    } else {

        sam_ref_ch.set{ indel_realigned_ch }
        
    }

    SAMTOOLS_SORT(indel_realigned_ch)
    | FILTER_BAM
    
    SAMTOOLS_FILTERED_INDEX(FILTER_BAM.out.bam_ch)
    | combine(ref)
    | SAMTOOLS_PILEUP
    | BCFTOOLS_CALL

    emit:
    BCFTOOLS_CALL.out.called_ch
}