nextflow.enable.dsl=2

// nextflow run scripts/main.nf --ref /home/ref.txt --program bwa --domapping True --human False --pairedend True --maxinsertsize 1000 --mininsertsize 50 --ssahaquality 30 --maprepeats False --GATK True --markdup True --detectOverlaps False --pseudosequence True --incref True --indels True --quality 50 --mapq 20 --depth 8 --stranddepth 3 --anomolous True --BAQ True --circular True --ratio 0.8 --prior 0.001 --call c --force False --filter 1 --tabfile False --alnfile False --raxml False --model GTRGAMMA --bootstrap 100 --keep False --LSF True --LSFQ normal --mem 5 --nodes 20 --dirty False --mapfiles "/home/ref.txt,/home/ref2.txt" -process.echo

include { mpfs } from './sub-workflows/makepileup_from_sam.nf'
include { MAPPING } from './modules/MAPPING.nf'

process LOG_COMMANDLINE {
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

process RANDOM_NAME {
    output:
    val tmpname

    exec:
    chars = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
    rnd = new Random()
    length = rnd.nextInt(3) + 8
    tmpname = 'tmp' + (1..length).collect { chars[random.nextInt(chars.length())] }.join('')
}

workflow {
    log_ch = LOG_COMMANDLINE()
    
    ziplist = [:].withDefault { [] }
    bamlist = [:]
    pools = []
    poolsort = []

    tmpname = RANDOM_NAME()

    (ziplist, bamlist, poolsort, pools) = MAPPING (Channel.from(params.mapfiles.split(',')).flatten(), tmpname, ziplist, bamlist, poolsort, pools)

    pools = pools.take(1)
    bamlist = bamlist.take(1)
    ziplist = ziplist.take(1)
    poolsort = poolsort.take(1)

    // pools.view()
    // bamlist.view()
    // poolsort.view()
    // ziplist.view()
}
