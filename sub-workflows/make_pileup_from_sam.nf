include { MARK_DUPLICATES; SEQUENCE_DICT                } from '../modules/picard.nf'
include { INDEX_REF
          SAMTOOLS_SORT_BAM_AND_MAKE_HEADER; 
          SAMTOOLS_SORT;
          SAMTOOLS_SORT as SAMTOOLS_SORT_RAW_MAPPING;
          SAMTOOLS_INDEX as SAMTOOLS_RAW_INDEX;
          SAMTOOLS_INDEX as SAMTOOLS_FILTERED_INDEX; 
          SAMTOOLS_MERGE;
          SAMTOOLS_PILEUP                               } from '../modules/samtools.nf'
include { INDEL_REALIGNMENT;
          ADD_READGROUP                                 } from '../modules/gatk_indel_realignment.nf'
include { FILTER_BAM                                    } from '../modules/filter_bam.nf'
include { BCFTOOLS_CALL                                 } from '../modules/bcftools.nf'

workflow MAKE_PILEUP_FROM_SAM {

    take:
    mapped_sam_ch
    ref

    main:

    if (params.markdup) {

        SAMTOOLS_SORT_RAW_MAPPING(mapped_sam_ch)
        | MARK_DUPLICATES

        MARK_DUPLICATES.out.deduped_ch
        | set { deduped_ch }

    } else {

        mapped_sam_ch.set{ deduped_ch }

    }

    SAMTOOLS_SORT_BAM_AND_MAKE_HEADER(deduped_ch)

    Map<String, String> program_label = [
        "BWA": "BWA MEM",
        "SSAHA": "SSAHA",
        "SMALT": "SMALT"
    ]

    ADD_READGROUP(
        SAMTOOLS_SORT_BAM_AND_MAKE_HEADER.out.header_ch,
        Channel.value(program_label[params.program])
    )
    | set { formatted_header_ch }
    
    SAMTOOLS_SORT_BAM_AND_MAKE_HEADER.out.bam_ch.join(formatted_header_ch)
    | SAMTOOLS_MERGE
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
    | set { sorted_indel_ch }

    FILTER_BAM(sorted_indel_ch)

    SAMTOOLS_FILTERED_INDEX(FILTER_BAM.out.bam_ch)
    | combine(ref)
    | SAMTOOLS_PILEUP
    | BCFTOOLS_CALL

    emit:
    BCFTOOLS_CALL.out.called_ch
}