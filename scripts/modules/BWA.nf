process BWA {
    input:
    path tmphead_sam
    val name

    output:
    path tmphead_sam

    script:
    """
    echo '@RG\\tID:${name}\\tCN:Sanger\\tDT:\$now\\tPG:BWA MEM\\tPL:ILLUMINA\\tSM:${name}' >> ${tmphead_sam}
    """
}