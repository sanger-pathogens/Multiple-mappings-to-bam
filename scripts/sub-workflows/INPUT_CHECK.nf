workflow INPUT_CHECK {
    take:
    tmpname

    main:
    ziplist = [:].withDefault { [] }
    poolsort = []
    dirs = []
    list_mp = params.mapfiles.split(',').collect()
    files = Channel.from(params.mapfiles.split(','))
        .map { row -> 
            process_inputs (row, tmpname, ziplist, poolsort, list_mp)
        }
        .map { pools, ziplist, bamlist, poolsort -> [pools, ziplist, bamlist, poolsort] }

    (zl, bl, ps) = PROCESS_LISTS(files.take(1))

    emit:
    files
    zl
    bl
    ps

}

process PROCESS_LISTS {
    input:
    val files

    output:
    val zl
    val bl
    val ps

    exec:
    zl = files[1]
    bl = files[2]
    ps = files[3]

    if (ps.size() == 0) {
        exit 1, "No valid input file"
    }
}

def process_inputs (String pool, String tmpname, Map<String, List<String>> ziplist, List<String> poolsort, List<String> list_mp) {
    bamlist = ""
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

            key = pool.replace("#", "_")
            bamlist=pool
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
            file1 = originalfastqdir+pool[0..-3]+"_2"+nonhumanpool+".fastq"
            file2 = originalfastqdir+pool[0..-3]+"_2"+nonhumanpool+".fastq.gz"
            if (!list_mp.contains(file1) && !list_mp.contains(file2)) {
                println "File "+pool+"_2.fastq not found! Treating "+pool+" as unpaired..."
                pairedend=false
            } else {
                pool = pool[0..-3]
            }
        } else if (pool[-2..-1] == '_2') {
            file1 = originalfastqdir+pool[0..-3]+"_1"+nonhumanpool+".fastq"
            file2 = originalfastqdir+pool[0..-3]+"_1"+nonhumanpool+".fastq.gz"
            if (!list_mp.contains(file1) && !list_mp.contains(file2)) {
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

    return [pools, ziplist, bamlist, poolsort]
}