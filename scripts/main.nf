nextflow.enable.dsl=2

// nextflow run scripts/main.nf --ref /home/darshil/LED/Nextflow/Streptococcus_agalactiae_NGBS128_GCF_001552035_1.fa --program BWA --domapping True --human False --pairedend True --maxinsertsize 1000 --mininsertsize 50 --ssahaquality 30 --maprepeats False --GATK True --markdup True --detectOverlaps False --pseudosequence True --incref True --indels True --quality 50 --mapq 20 --depth 8 --stranddepth 3 --anomolous True --BAQ True --circular True --ratio 0.8 --prior 0.001 --call c --output Streptococcus_agalactiae_NGBS128_GCF_001552035_1_bwa --force False --filter 1 --tabfile False --alnfile False --raxml False --model GTRGAMMA --bootstrap 100 --keep False --LSF True --LSFQ normal --mem 5 --nodes 20 --dirty False --mapfiles "/home/darshil/LED/Nextflow/31663_7#10_1.fastq.gz,/home/darshil/LED/Nextflow/31663_7#10_2.fastq.gz,/home/darshil/LED/Nextflow/31663_7#12_1.fastq.gz,/home/darshil/LED/Nextflow/31663_7#12_2.fastq.gz" -process.echo

//MODULES
include { BWA_INDEX } from './modules/BWA.nf'
include { SMALT_INDEX } from './modules/SMALT.nf'
include { HANDLE_SEQUENCES } from './modules/HANDLE_SEQUENCES'

//SUBWORKFLOWS
include { INPUT_CHECK } from './sub-workflows/INPUT_CHECK.nf'
include { CALL_MAPPING } from './sub-workflows/CALL_MAPPING.nf'
include { MAKE_PILEUP_FROM_SAM } from './sub-workflows/MAKE_PILEUP_FROM_SAM.nf'
include { PSEUDOSEQUENCE_GENERATION } from './sub-workflows/PSEUDOSEQUENCE_GENERATION.nf'


process LOG_COMMANDLINE {
    label "cpu_1"
    label "mem_16"
    label "time_1"
    
    publishDir "${params.outdir}", mode: 'copy'
    output:
    path 'MM_command_*.txt'

    script:
    """
    timestamp=\$(date +"%Y-%m-%d.%H.%M.%S")

    output_file="MM_command_\$timestamp.txt"

    echo "nextflow run scripts/main.nf \
    --ref ${params.ref} \
    --program ${params.program} \
    --domapping ${params.domapping} \
    --human ${params.human} \
    --pairedend ${params.pairedend} \
    --maxinsertsize ${params.maxinsertsize} \
    --mininsertsize ${params.mininsertsize} \
    --ssahaquality ${params.ssahaquality} \
    --maprepeats ${params.maprepeats} \
    --GATK ${params.GATK} \
    --markdup ${params.markdup} \
    --detectOverlaps ${params.detectOverlaps} \
    --pseudosequence ${params.pseudosequence} \
    --incref ${params.incref} \
    --indels ${params.indels} \
    --quality ${params.quality} \
    --mapq ${params.mapq} \
    --depth ${params.depth} \
    --stranddepth ${params.stranddepth} \
    --anomolous ${params.anomolous} \
    --BAQ ${params.BAQ} \
    --circular ${params.circular} \
    --ratio ${params.ratio} \
    --prior ${params.prior} \
    --call ${params.call} \
    --force ${params.force} \
    --filter ${params.filter} \
    --tabfile ${params.tabfile} \
    --alnfile ${params.alnfile} \
    --raxml ${params.raxml} \
    --model ${params.model} \
    --bootstrap ${params.bootstrap} \
    --keep ${params.keep} \
    --LSF ${params.LSF} \
    --LSFQ ${params.LSFQ} \
    --mem ${params.mem} \
    --nodes ${params.nodes} \
    --dirty ${params.dirty} \
    --mapfiles ${params.mapfiles} \
    -process.echo" >> \$output_file
    """
}

workflow {
    log_ch = LOG_COMMANDLINE()

    ref = Channel.fromPath(params.ref)

    read_ch = INPUT_CHECK()

    (ref, ref_aln) = HANDLE_SEQUENCES(ref)
    ref = ref.collect()

    if (params.program == "BWA") {
        (index_ch, fai) = BWA_INDEX(ref)
    } else if (params.program == "SMALT") {
        (index_ch, fai) = SMALT_INDEX(ref)
    }
    
    CALL_MAPPING(read_ch, ref, fai, index_ch)
    | MAKE_PILEUP_FROM_SAM
    | set { called_ch }

    if (params.pseudosequence == true) {
        PSEUDOSEQUENCE_GENERATION(ref, called_ch)
    }

}