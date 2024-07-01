process MAPPING {
    input:
    val pool

    script:
    def originalfastqdir=""
    def nonhumanpool=""
    """
    pool=${pool}
    if [ ! -f "${pool}" ]; then
        echo "File ${pool} not found! Skipping..."
        exit 0
    fi

    file_extension="\${pool##*.}"

    if [[ "\${file_extension}" != "fastq" && "\${file_extension}" != "bam" && "\${file_extension}" != "gz" ]]; then
        echo "WARNING: Input file name is not .fastq, .bam, or .gz!"
    fi

    if [ "${pool: -1}" == "/" ]; then
        pool="${pool%/}"
    fi

    originalfastqdir=$(dirname "\$pool")
    """

}