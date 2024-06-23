nextflow.enable.dsl=2

// nextflow run scripts/main.nf --ref /home/ref.txt --program bwa --domapping True --human False --pairedend True --maxinsertsize 1000 --mininsertsize 50 --ssahaquality 30 --maprepeats False --GATK True --markdup True --detectOverlaps False --pseudosequence True --incref True --indels True --quality 50 --mapq 20 --depth 8 --stranddepth 3 --anomolous True --BAQ True --circular True --ratio 0.8 --prior 0.001 --call c --force False --filter 1 --tabfile False --alnfile False --raxml False --model GTRGAMMA --bootstrap 100 --keep False --LSF True --LSFQ normal --mem 5 --nodes 20 --dirty False --mapfiles "/home/ref.txt,/home/ref2.txt" -process.echo

include { runSsaha } from './modules/runssaha.nf'
include { runSMALT } from './modules/runsmalt.nf'
include { makePileupFromSam } from './modules/make_pileup_from_sam.nf'

process processArgs {

    shell:
    """
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
    echo "Mapfiles: ${params.mapfiles}"
    """
}

workflow {
    
}
