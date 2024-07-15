workflow INPUT_CHECK {
    take:
    tmpname

    main:
    ziplist = [:].withDefault { [] }
    bamlist = [:]
    poolsort = []
    // pools = []
    dirs = []

    files = Channel.from(params.mapfiles.split(','))
        .map { row -> 
            mapping (row, tmpname, ziplist, bamlist, poolsort, dirs)
        }
        .map { pools, ziplist, bamlist, poolsort, dirs -> [pools, ziplist, bamlist, poolsort, dirs] }

    
    make_dirs = PROCESS_DIRS(files.take(1))
    GENERATE_DIR(make_dirs.flatten())

}
process PROCESS_DIRS {
    input:
    val files

    output:
    val make_dirs

    exec:
    make_dirs = files[4]
}
process GENERATE_DIR {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val dir

    output:
    path "${dir}"

    script:
    """
    mkdir -p ${dir}
    """
}

def mapping (String pool, String tmpname, Map<String, String> ziplist, Map<String, String> bamlist, List<String> poolsort, List<String> dirs) {
    // println tmpname
    if (pool[-1] == '/') {
        pool = pool[0..-2]
    }
    println pool
    filetype = '.' + pool.split('\\.')[-1]
    println "filetype ${filetype}"

    if (!['.fastq', '.bam', '.gz'].contains(filetype)) {
        println "WARNING: Input file name is not .fastq or .bam!"
    }

    bam = ""
    originalfastqdir = ''
    if (pool.split('/').size() > 1) {
        originalfastqdir = pool.split('/').take(pool.split('/').size() - 1).join('/') + '/'
    }

    nonhumanpool = ''
    if (pool.contains("_nonhuman")) {
        nonhumanpool = "_nonhuman"
        pool = pool.replace("_nonhuman", "")
    }

    is_zipped = false
    unbammed = false
    unzipped = false
    if (params.program != 'bwa' && pool.split('\\.')[-1] == "gz" && pool.split('\\.')[-2] == "fastq") {
        // command = "mkdir -p ${tmpname}_unzipped"
        // process = command.execute()
        // process.waitFor()
        dirs << "${tmpname}_unzipped"

        if (params.pairedend) {
            if (pool.split('\\.')[-3].endsWith("_1") || pool.split('\\.')[-3].endsWith("_2")) {
                ziplist[(pool.split('/')[-1].split('\\.')[0..-3]).join('.')[0..-3]] << pool
            } else {
                ziplist[(pool.split('/')[-1].split('\\.')[0..-3]).join('.')] << pool
            }
        } else {
            ziplist[(pool.split('/')[-1].split('\\.')[0..-3]).join('.')] << pool
        }
        
        pool = "${tmpname}_unzipped/" + (pool.split('/')[-1].split('\\.')[0..-2]).join('.')
        is_zipped = true
    } else if (params.program == 'bwa' && pool.split('\\.')[-1] == "gz") {
        is_zipped = true
        pool = originalfastqdir + (pool.split('/')[-1].split('\\.')[0..-2]).join('.')
    } else if (pool.split('\\.')[-1] == "bam") {
        if (params.domapping) {
            // command = "mkdir -p ${tmpname}_unbammed"
            // process = command.execute()
            // process.waitFor()
            dirs << "${tmpname}_unbammed"

            key = pool.replace("#", "_")
            bamlist[key.split('/')[-1].split('\\.')[0..-2].join('.')]=pool
            pool = pool.replace("#", "_")
            pool = "${tmpname}_unbammed/" + pool.split('/')[-1].split('\\.')[0..-2].join('.') + ".bam"
        } else {
            bam = pool
        }
    }

    fastqdir=''

    if (pool[-1] == '/') {
        pool = pool[0..-2]
    }
    if (pool.split('/').size() > 1) {
        fastqdir = pool.split('/')[0..-2].join('/') + "/"
    }

    pool = pool.split('/')[-1].split('\\.')[0..-2].join('.')
    pairedend = ""
    if (params.pairedend && filetype != ".bam") {
        if (pool[-2..-1] == '_1') {
            file1 = new File(originalfastqdir+pool[0..-3]+"_2"+nonhumanpool+".fastq")
            file2 = new File(originalfastqdir+pool[0..-3]+"_2"+nonhumanpool+".fastq.gz")
            if (!file1.exists() && !file2.exists()) {
                println "File "+pool+"_2.fastq not found! Treating "+pool+" as unpaired..."
                pairedend=false
            } else {
                pool = pool[0..-3]
            }
        } else if (pool[-2..-1] == '_2') {
            file1 = new File(originalfastqdir+pool[0..-3]+"_1"+nonhumanpool+".fastq")
            file2 = new File(originalfastqdir+pool[0..-3]+"_1"+nonhumanpool+".fastq.gz")
            if (!file1.exists() && !file2.exists()) {
                println "File "+pool+"_1.fastq not found! Treating "+pool+" as unpaired..."
                pairedend=false
            } else {
                pool = pool[0..-3]
            }
        } else {
            println "Not a typical paired-end name format! Treating "+pool+" as unpaired..."
            pairedend=false
        }
    }

    name = pool

    if (params.diroutput=="") {
        pool = pool + "_" + params.program
    } else {
        pool = pool + "_" + params.diroutput
    }

    if (pool in poolsort) {
        return
    }
    println pool+'...'

    // command = "mkdir -p ${pool}"
    // process = command.execute()
    // process.waitFor()
    dirs << "${pool}"

    pools = [:]
    pools["runname"] = pool
    pools["name"] = name
    pools["fastqdir"] = fastqdir
    pools["filetype"] = filetype
    pools["pairedend"] = pairedend
    pools["is_zipped"] = is_zipped

    if (bam != "") {
        pools["domapping"] = false
        pools["bam"] = bam
    } else {
        pools["domapping"] = true
    }
    // println "pools "+ pools
    // pools << pools
    // println "pools "+ pools
    poolsort << pool

    println "ok"

    return [pools, ziplist, bamlist, poolsort, dirs]
}