def log_commandline() {
    def timestamp = new Date().format("yyyy-MM-dd.HH.mm.ss")
    def output_file = "MM_command_${timestamp}.txt"

    log.info("""
Date: ${timestamp}

Nextflow Command:
nextflow run scripts/main.nf \\
    --ref ${params.ref} \\
    --program ${params.program} \\
    --domapping ${params.domapping} \\
    --human ${params.human} \\
    --pairedend ${params.pairedend} \\
    --maxinsertsize ${params.maxinsertsize} \\
    --mininsertsize ${params.mininsertsize} \\
    --ssahaquality ${params.ssahaquality} \\
    --maprepeats ${params.maprepeats} \\
    --GATK ${params.GATK} \\
    --markdup ${params.markdup} \\
    --detectOverlaps ${params.detectOverlaps} \\
    --pseudosequence ${params.pseudosequence} \\
    --incref ${params.incref} \\
    --indels ${params.indels} \\
    --quality ${params.quality} \\
    --mapq ${params.mapq} \\
    --depth ${params.depth} \\
    --stranddepth ${params.stranddepth} \\
    --dontuseanomolous ${params.dontuseanomolous} \\
    --BAQ ${params.BAQ} \\
    --circular ${params.circular} \\
    --ratio ${params.ratio} \\
    --prior ${params.prior} \\
    --call ${params.call} \\
    --force ${params.force} \\
    --filter ${params.filter} \\
    --tabfile ${params.tabfile} \\
    --alnfile ${params.alnfile} \\
    --raxml ${params.raxml} \\
    --model ${params.model} \\
    --bootstrap ${params.bootstrap} \\
    --keep ${params.keep} \\
    --LSF ${params.LSF} \\
    --LSFQ ${params.LSFQ} \\
    --mem ${params.mem} \\
    --nodes ${params.nodes} \\
    --dirty ${params.dirty} \\
    --read_dir ${params.read_dir}
""")
}