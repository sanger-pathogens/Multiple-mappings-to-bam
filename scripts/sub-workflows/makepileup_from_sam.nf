include { SORT_AND_MARK_DUPLICATES } from './../modules/SORT_AND_MARK_DUPLICATES.nf'
include { SAMTOOLS_SORT;SAMTOOLS_SORT_1;SAMTOOLS_INDEX;SAMTOOLS_MERGE } from './../modules/SAMTOOLS_SORT.nf'
include { SMALT } from './../modules/SMALT.nf'
include { BWA } from './../modules/BWA.nf'
include { INDEL_REALIGNMENT } from './../modules/INDEL_REALIGNMENT.nf'
include { FILTER_BAM } from './../modules/FILTER_BAM.nf'
include { PILEUP } from './../modules/PILEUP.nf'
include { BCFTOOLS_CALL } from './../modules/BCFTOOLS_CALL.nf'

import java.math.BigDecimal

workflow mpfs {

    take:
    name_tmp1_bam      //tuple(name, path tmp1_bam)
    cmdline

    main:

    ref = channel.fromPath(params.ref)

    if (params.markdup == "true") {
        (tmp1_bam, matrix_file_ch) = SORT_AND_MARK_DUPLICATES(name_tmp1_bam)
    }

    (name_bam, tmphead_sam, bam_bai) = SAMTOOLS_SORT(name_tmp1_bam)

    if (params.program == "smalt" || params.program == "SMALT") {
        tmphead_sam = SMALT(tmphead_sam, cmdline)
    } else if (params.program == "bwa" || params.program == "BWA") {
        tmphead_sam = BWA(tmphead_sam, name)
    }

    (tmp1_bam, tmphead_bam) = SAMTOOLS_MERGE(name_bam, tmphead_sam[1])

    if (params.GATK == 'true') {
        tmp1_bam = INDEL_REALIGNMENT(tmp1_bam, ref)
    }

    tmp_bam = SAMTOOLS_SORT_1(tmp1_bam)

    (name_bam, filter_bam_ch) = FILTER_BAM(tmp_bam)

    bam_bai = SAMTOOLS_INDEX(name_bam[1])

    tmp_mpileup = PILEUP(ref, name_bam)

    (name_ploidy, name_bcf, name_variant_bcf, name_bcf_csi, name_variant_bcf_csi) = BCFTOOLS_CALL(tmp_mpileup)

    if (params.pseudosequence == "true") {
        pseudosequence = PSEUDOSEQUENCE(name_bcf, name_bam[1])
    }

    emit:
    // All the files generated above
    matrix_file_ch
    name_bam
    filter_bam_ch
    bam_bai
    tmp_mpileup
    name_ploidy
    name_bcf
    name_variant_bcf
    name_variant_bcf_csi
    tmp1_bam
    tmphead_sam
    tmphead_bam
    tmp_bam
    pseudosequence

}