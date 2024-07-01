process SMALT {
    input:
    path tmphead_sam        //${runname}/tmphead.sam
    val cmdline

    output:
    path tmphead_sam

    script:
    """
    now=\$(date +'%Y-%m-%dT%H:%M:%S')
    echo "@RG\\tID:${name}\\tCN:Sanger\\tDT:\$now\\tPG:SMALT\\tPL:ILLUMINA\\tSM:${name}" >> ${tmphead_sam}
    if [ ${params.domapping} ] && [ ${newsmalt} = "false" ]
    then
        smaltversion=\$( smalt version | grep Version | awk '{print \$2}' )
        echo "@PG\\tID:SMALT\\tPN:SMALT\\tCL:${cmdline}\\tVN:\$smaltversion" >> ${tmphead_sam}
    fi
    """
}