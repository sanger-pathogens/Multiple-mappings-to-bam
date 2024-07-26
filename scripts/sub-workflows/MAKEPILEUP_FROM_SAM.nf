include { SORT_AND_MARK_DUPLICATES } from './../modules/SORT_AND_MARK_DUPLICATES.nf'
include { SAMTOOLS_SORT;SAMTOOLS_SORT_1;SAMTOOLS_INDEX;SAMTOOLS_MERGE } from './../modules/SAMTOOLS_SORT.nf'
include { SMALT } from './../modules/SMALT.nf'
include { BWA } from './../modules/BWA.nf'
include { INDEL_REALIGNMENT } from './../modules/INDEL_REALIGNMENT.nf'
include { FILTER_BAM } from './../modules/FILTER_BAM.nf'
include { PILEUP } from './../modules/PILEUP.nf'
include { BCFTOOLS_CALL } from './../modules/BCFTOOLS_CALL.nf'
include { PSEUDOSEQUENCE } from './../modules/PSEUDOSEQUENCE.nf'

import java.math.BigDecimal

workflow MAKEPILEUP_FROM_SAM {

    take:
    files

    main:

    ref = channel.fromPath(params.ref)

    if (params.markdup == "true") {
        (files, matrix_file_ch) = SORT_AND_MARK_DUPLICATES(files)
    }

    files = SAMTOOLS_SORT(files)

    if (params.program == "SMALT") {
        files = SMALT(files)
    } else if (params.program == "BWA") {
        files = BWA(files)
    }

    (files, tmphead_bam) = SAMTOOLS_MERGE(files)

    if (params.GATK == true) {
        files = INDEL_REALIGNMENT(files, ref)
    }

    files = SAMTOOLS_SORT_1(files)

    (files, filter_bam_ch) = FILTER_BAM(files)

    files = SAMTOOLS_INDEX(files)

    files = PILEUP(files, ref)

    (files, name_ploidy, name_variant_bcf, name_bcf_csi, name_variant_bcf_csi) = BCFTOOLS_CALL(files)

    if (params.pseudosequence == true) {
        pseudosequence = PSEUDOSEQUENCE(files)
    }

    emit:
    // All the files generated above
    files
    // tmphead_bam
    // matrix_file_ch
    // pseudosequence
    // filter_bam_ch
    // name_ploidy
    // name_variant_bcf
    // name_bcf_csi
    // name_variant_bcf_csi
}