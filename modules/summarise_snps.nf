process SUMMARISE_SNPS {
    label "cpu_1"
    label "mem_1"
    label "time_1"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "*.out"

    container 'quay.io/ssd28/gsoc-experimental/summarise_snps:0.0.3'

    input:
    tuple path(output_aln), path(ref)

    output:
    tuple path("${output_aln}.out"), path("${output_aln}_summary.out")

    script:
    summarystring = "summarise_snps.py -g -w -r "+ "${ref.baseName}" + " -o " + "${output_aln}" + " -i " + output_aln

    if (params.embl != "") {
        summarystring = summarystring + " -e "+ params.embl
    }
    if (params.alnfile == true) {
        summarystring = summarystring + " -a"
    }
    if (params.tabfile == true) {
        summarystring = summarystring + " -t"
    }
    if (params.raxml == true) {
        summarystring = summarystring + " -p -l -b ${params.bootstrap}"
    }

    """
    ${summarystring}
    """
}