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
    files

    main:

    ref = channel.fromPath(params.ref)

    if (params.markdup == "true") {
        (files, matrix_file_ch) = SORT_AND_MARK_DUPLICATES(files)
    }

    files = SAMTOOLS_SORT_BAM_AND_MAKE_HEADER(files)

    if (params.program == "SMALT") {
        files = FORMAT_SMALT_HEADER(files)
    } else if (params.program == "BWA") {
        files = FORMAT_BWA_HEADER(files)
    }

    (files, tmphead_bam) = SAMTOOLS_MERGE(files)
    files_ref = files.combine(ref)
    if (params.GATK == true) {
        files = INDEL_REALIGNMENT(files_ref)
    }

    files = SAMTOOLS_SORT(files)

    (files, filter_bam_ch) = FILTER_BAM(files)

    files = SAMTOOLS_INDEX(files)
    files_ref = files.combine(ref)
    files = PILEUP(files_ref)

    (called_ch, name_ploidy, name_variant_bcf, name_bcf_csi, name_variant_bcf_csi) = BCFTOOLS_CALL(files)

    emit:
    called_ch
}