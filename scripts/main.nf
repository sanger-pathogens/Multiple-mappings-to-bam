nextflow.enable.dsl=2

include { validateInputs } from './validate_input.nf'

params.mapfiles = ""
params.ref = ""
params.program = "BWA"
params.domapping = true
params.human = false
params.pairedend = true
params.maxinsertsize = 1000
params.mininsertsize = 50
params.ssahaquality = 30
params.maprepeats = false
params.GATK = true
params.markdup = true
params.detectOverlaps = false
params.pseudosequence = true
params.incref = true
params.indels = true
params.quality = 50
params.mapq = 20
params.depth = 8
params.stranddepth = 3
params.anomolous = true
params.BAQ = true
params.circular = true
params.ratio = 0.8
params.prior = 0.001
params.call = "c"
params.embl = ""
params.output = " "
params.diroutput = ""
params.force = false
params.filter = "1"
params.tabfile = false
params.alnfile = false
params.raxml = false
params.model = "GTRGAMMA"
params.bootstrap = 100
params.keep = false
params.LSF = true
params.LSFQ = "normal"
params.mem = 5
params.nodes = 20
params.dirty = false

process processArgs {

    shell:
    """
    echo "Mapfiles: ${params.mapfiles}"
    echo "Ref: ${params.ref}"
    echo "Program: ${params.program}"
    echo "Domapping: ${params.domapping}"
    echo "Human: ${params.human}"
    echo "Pairedend: ${params.pairedend}"
    echo "Maxinsertsize: ${params.maxinsertsize}"
    echo "Mininsertsize: ${params.mininsertsize}"
    echo "Ssahaquality: ${params.ssahaquality}"
    echo "Maprepeats: ${params.maprepeats}"
    echo "GATK: ${params.GATK}"
    echo "Markdup: ${params.markdup}"
    echo "DetectOverlaps: ${params.detectOverlaps}"
    echo "Pseudosequence: ${params.pseudosequence}"
    echo "Incref: ${params.incref}"
    echo "Indels: ${params.indels}"
    echo "Quality: ${params.quality}"
    echo "Mapq: ${params.mapq}"
    echo "Depth: ${params.depth}"
    echo "Stranddepth: ${params.stranddepth}"
    echo "Anomolous: ${params.anomolous}"
    echo "BAQ: ${params.BAQ}"
    echo "Circular: ${params.circular}"
    echo "Ratio: ${params.ratio}"
    echo "Prior: ${params.prior}"
    echo "Call: ${params.call}"
    echo "Force: ${params.force}"
    echo "Filter: ${params.filter}"
    echo "Tabfile: ${params.tabfile}"
    echo "Alnfile: ${params.alnfile}"
    echo "Raxml: ${params.raxml}"
    echo "Model: ${params.model}"
    echo "Bootstrap: ${params.bootstrap}"
    echo "Keep: ${params.keep}"
    echo "LSF: ${params.LSF}"
    echo "LSFQ: ${params.LSFQ}"
    echo "Mem: ${params.mem}"
    echo "Nodes: ${params.nodes}"
    echo "Dirty: ${params.dirty}"
    """
}

workflow {
    validateInputs()
    processArgs()
}