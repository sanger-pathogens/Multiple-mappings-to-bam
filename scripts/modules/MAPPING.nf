process MAPPING {
    publishDir "results", mode: 'copy'

    input:
    val pool
    val tmpname
    val ziplist
    val bamlist
    val poolsort
    val pools

    output:
    val ziplist
    val bamlist
    val poolsort
    val pools
    path "${pool}"
    path "${tmpname}_unbammed", optional: true
    path "${tmpname}_unzipped", optional: true

    script:
    println pool
    if (pool[-1] == '/') {
        pool = pool[0..-2]
    }
    filetype = '.' + pool.split('\\.')[-1]

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
        unzipped = true

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
            unbammed = true

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
        println "returning"
        println "pool " + pool
        println "poolsort " + poolsort
        return
    }
    println pool+'...'

    // command = "mkdir -p ${pool}"
    // process = command.execute()
    // process.waitFor()

    mp = [:]
    mp["runname"] = pool
    mp["name"] = name
    mp["fastqdir"] = fastqdir
    mp["filetype"] = filetype
    mp["pairedend"] = pairedend
    mp["is_zipped"] = is_zipped

    if (bam != "") {
        mp["domapping"] = false
        mp["bam"] = bam
    } else {
        mp["domapping"] = true
    }
    // println "pools "+ pools
    pools << mp
    // println "pools "+ pools
    poolsort << pool

    println "ok"
    """
    if [ "${unbammed}" = "true" ]; then
        mkdir -p "${tmpname}_unbammed"
    elif [ "${unzipped}" = "true" ]; then
        mkdir -p "${tmpname}_unzipped"
    fi
    mkdir -p "${pool}"
    """
}