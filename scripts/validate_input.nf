process validateInputs {
    script:
    """
    echo "Validating input parameters"
    if [ -z "${params.ref}" ]; then
        echo "No reference dna file (-r) selected!"
        exit 1
    fi

    if [ ${params.maxinsertsize} -gt 100000 ] || [ ${params.maxinsertsize} -lt 10 ]; then
        echo "Maximum insert size (-i) must be between 10 and 100,000"
        exit 1
    fi

    if [ ${params.maxinsertsize} -lt ${params.mininsertsize} ]; then
        echo "Minimum insert size (-j) must be smaller than maximum insert size (-i). Currently -i=${params.maxinsertsize} and -j=${params.mininsertsize}"
        exit 1
    fi

    if [ ${params.mininsertsize} -gt 10000 ] || [ ${params.mininsertsize} -lt 10 ]; then
        echo "Minimum insert size (-j) must be between 10 and 10,000"
        exit 1
    fi

    if [ ${params.quality} -gt 99 ] || [ ${params.quality} -lt 1 ]; then
        echo "Ssaha mapping quality score (-q) must be between 1 and 100"
        exit 1
    fi

    if [ "${params.program}" = "BWA" ] && ( [ ${params.mapq} -gt 30 ] || [ ${params.mapq} -lt 0 ] ); then
        echo "Mapping quality score (-Q) must be between 0 and 30 for bwa"
        exit 1
    fi

    if [ "${params.program}" != "BWA" ] && ( [ ${params.mapq} -gt 60 ] || [ ${params.mapq} -lt 0 ] ); then
        echo "Mapping quality score (-Q) must be between 0 and 60 for ${params.program}"
        exit 1
    fi

    if [ ${params.ratio} -gt 1 ] || [ ${params.ratio} -lt 0 ]; then
        echo "SNP/site mapping quality ratio cutoff (-R) must be between 0 and 1"
        exit 1
    fi

    if [ ${params.prior} -gt 1 ] || [ ${params.prior} -lt 0 ]; then
        echo "Estimated mutation rate (-P) must be between 0 and 1"
        exit 1
    fi

    if [ ${params.mem} -gt 30 ] || [ ${params.mem} -lt 0 ]; then
        echo "Memory requirement (-M) must be between 0 and 30Gb"
        exit 1
    fi



    if [ "${params.force}" = "false" ] && [ -f "${params.output}.aln" ]; then
        echo "Output files with chosen prefix already exist. Use force option to overwrite."
        exit 1
    fi
    """
}