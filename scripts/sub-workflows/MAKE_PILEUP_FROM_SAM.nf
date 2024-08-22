include { SORT_AND_MARK_DUPLICATES } from './../modules/SORT_AND_MARK_DUPLICATES.nf'
include { SAMTOOLS_SORT_BAM_AND_MAKE_HEADER; SAMTOOLS_SORT; SAMTOOLS_INDEX; SAMTOOLS_MERGE } from './../modules/SAMTOOLS_SORT.nf'
include { FORMAT_SMALT_HEADER } from './../modules/SMALT.nf'
include { FORMAT_BWA_HEADER } from './../modules/BWA.nf'
include { INDEL_REALIGNMENT } from './../modules/INDEL_REALIGNMENT.nf'
include { FILTER_BAM } from './../modules/FILTER_BAM.nf'
include { PILEUP } from './../modules/PILEUP.nf'
include { BCFTOOLS_CALL } from './../modules/BCFTOOLS_CALL.nf'

import java.math.BigDecimal

workflow MAKE_PILEUP_FROM_SAM {

    take:
    mapped_sam_ch

    main:

    ref = channel.fromPath(params.ref)

    if (params.markdup == "true") {
        (deduped_ch, matrix_file_ch) = SORT_AND_MARK_DUPLICATES(mapped_sam_ch)
    } else {
        mapped_sam_ch.set{ deduped_ch }
    }

    sorted_header_ch = SAMTOOLS_SORT_BAM_AND_MAKE_HEADER(deduped_ch)

    if (params.program == "SMALT") {
        formatted_sam_ch = FORMAT_SMALT_HEADER(sorted_header_ch)
    } else if (params.program == "BWA") {
        formatted_sam_ch = FORMAT_BWA_HEADER(sorted_header_ch)
    }
    

    (merged_ch, tmphead_bam) = SAMTOOLS_MERGE(formatted_sam_ch)


    sam_ref_ch = merged_ch.combine(ref)


    if (params.GATK == true) {
        indel_realigned_ch = INDEL_REALIGNMENT(sam_ref_ch)
    } else {
        sam_ref_ch.set{ indel_realigned_ch }
    }

    sorted_indel_ch = SAMTOOLS_SORT(indel_realigned_ch)

    (filtered_ch, filter_bam_ch) = FILTER_BAM(sorted_indel_ch)

    indexed_ch = SAMTOOLS_INDEX(filtered_ch)
    indexed_plus_ref_ch = indexed_ch.combine(ref)
    pileup_ch = PILEUP(indexed_plus_ref_ch)

    (called_ch, name_ploidy, name_variant_bcf, name_bcf_csi, name_variant_bcf_csi) = BCFTOOLS_CALL(pileup_ch)

    emit:
    called_ch
}