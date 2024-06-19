def getSamtoolsVersion() {
    try {
        def samtools_version_text = "samtools --version".execute().text
        def samtools_version_firstline = samtools_version_text.readLines().first()
        def version = samtools_version_firstline.split()[1]
        return version
    } catch (Exception e) {
        println e
        return 0
    }
}

process makePileupFromSam {
    input:
    path bashfile
    val runname
    val name

    script:
    def version = getSamtoolsVersion()

    """
    chevron=""
    suffix=""
    anomolous=""
    if [ ${version} -gt 1.2 ]
    then
        chevron=">"
        suffix=".bam"
    fi

    if [ ${params.markdup} ]
    then
        echo "samtools sort ${runname}/tmp1.bam ${chevron} ${runname}/tmpsort${suffix}" >> ${bashfile}
        echo "picard MarkDuplicates INPUT=${runname}/tmpsort.bam OUTPUT=${runname}/tmp1.bam METRICS_FILE=${runname}/${name}_metrics.txt" >> ${bashfile}
        echo "rm ${runname}/tmpsort.bam" >> ${bashfile}
    fi

    echo "samtools sort ${runname}/tmp1.bam ${chevron} ${runname}/${name}${suffix}" >> ${bashfile}
    echo "samtools index ${runname}/${name}.bam" >> ${bashfile}
    echo "rm ${runname}/tmp1.bam" >> ${bashfile}

    echo "samtools view -H ${runname}/${name}.bam | sed 's/SO:unknown/SO:coordinate/g' | sed 's/\\\\x00//g' > ${runname}/tmphead.sam" >> ${bashfile}

    now=$(date +'%Y-%m-%dT%H:%M:%S')

    if [ "${params.program}" = "smalt" ] || [ "${params.program}" = "SMALT" ]
    then
        echo "echo \"@RG\\tID:${name}\\tCN:Sanger\\tDT:\$now\\tPG:SMALT\\tPL:ILLUMINA\\tSM:${name}\" >> ${runname}/tmphead.sam" >> ${bashfile}
        if [ ${params.domapping} ] && [ ! ${params.newsmalt} ]
        then
            echo "smaltversion=$( smalt version | grep Version | awk '{print \$2}' )" >> ${bashfile}
            echo "echo \"@PG\\tID:SMALT\\tPN:SMALT\\tCL:${params.cmdline}\\tVN:\$smaltversion\" >> ${runname}/tmphead.sam" >> ${bashfile}
        fi
    elif [ "${params.program}" = "bwa" ] || [ "${params.program}" = "BWA" ]
    then
        echo "echo '@RG\\tID:${name}\\tCN:Sanger\\tDT:\$now\\tPG:BWA MEM\\tPL:ILLUMINA\\tSM:${name}' >> ${runname}/tmphead.sam" >> ${bashfile}
    fi

    echo "samtools view -b -o ${runname}/tmphead.bam -H ${runname}/${name}.bam" >> ${bashfile}
    echo "samtools merge -c -p -f -r -h ${runname}/tmphead.sam ${runname}/tmp.bam ${runname}/${name}.bam ${runname}/tmphead.bam" >> ${bashfile}
    echo "mv ${runname}/tmp.bam ${runname}/tmp1.bam" >> ${bashfile}
    echo "rm ${runname}/${name}.bam" >> ${bashfile}

    if [ ${params.GATK} ]
    then
        echo "samtools index ${runname}/tmp1.bam" >> ${bashfile}
        echo "cp ${params.ref} ${runname}/tmpref.fa" >> ${bashfile}
        echo "samtools faidx ${runname}/tmpref.fa" >> ${bashfile}
        echo "picard CreateSequenceDictionary R=${runname}/tmpref.fa O=${runname}/tmpref.dict" >> ${bashfile}
        echo "gatk -I ${runname}/tmp1.bam -R ${runname}/tmpref.fa -T RealignerTargetCreator -o ${runname}/tmp.intervals" >> ${bashfile}
        echo "gatk -I ${runname}/tmp1.bam -R ${runname}/tmpref.fa -T IndelRealigner --filter_bases_not_stored -targetIntervals ${runname}/tmp.intervals -o ${runname}/tmp.bam" >> ${bashfile}
        echo "mv ${runname}/tmp.bam ${runname}/tmp1.bam" >> ${bashfile}
        echo "rm ${runname}/tmp1.bam.bai ${runname}/tmpref.* ${runname}/tmp.intervals ${runname}/tmphead.*" >> ${bashfile}
    fi

    echo "samtools sort ${runname}/tmp1.bam ${chevron} ${runname}/tmp${suffix}" >> ${bashfile}
    echo "rm ${runname}/tmp1.bam" >> ${bashfile}

    if [ "${params.filter}" = "1" ]
    then
        echo "mv ${runname}/tmp.bam ${runname}/${name}.bam" >> ${bashfile}
    elif [ "${params.filter}" = "2" ]
    then
        echo "samtools view -F 4 -b -o ${runname}/${name}.bam ${runname}/tmp.bam"
    elif [ "${params.filter}" = "3" ]
    then
        echo "samtools view -f 2 -b -o ${runname}/${name}.bam ${runname}/tmp.bam"
    elif [ "${params.filter}" = "4" ]
    then
        echo "samtools view -F 4 -b -o ${runname}/${name}.bam ${runname}/tmp.bam && samtools view -f 4 -b -o ${runname}/${name}_unmapped.bam ${runname}/tmp.bam"
    elif [ "${params.filter}" = "5" ]
    then
        echo "samtools view -f 2 -b -o ${runname}/${name}.bam ${runname}/tmp.bam && samtools view -F 2 -b -o ${runname}/${name}_unpaired.bam ${runname}/tmp.bam"
    fi

    echo "samtools index ${runname}/${name}.bam" >> ${bashfile}

    if [ ${params.anomolous} ]
    then
        anomolous=" -A "
    else
        anomolous=""
    fi

    echo "echo \"${name}    1 > ${runname}/${name}.ploidy\"" >> ${bashfile}

    BAQ=""
    overlaps=""

    if [ "${params.program}" = "bwa" ] || [ "${params.program}" = "BWA" ]
    then
        if [ ${params.BAQ} = "false" ]
        then
            BAQ="-B"
        else
            BAQ=""
        fi
        if [ ${params.detectOverlaps} ]
        then
            overlaps=""
        else
            overlaps="-x"
        fi
        echo "samtools mpileup -t DP,DP4 -C 50 -L 1000 -d 1000 -m ${params.depth} \$anomolous \$BAQ \$overlaps -ugf ${params.ref} ${runname}/${name}.bam > ${runname}/tmp.mpileup" >> ${bashfile}
    else
        if [ ${params.BAQ} = "true" ]
        then
            BAQ=""
        else
            BAQ="-B"
        fi
        if [ ${params.detectOverlaps} = "true" ]
        then
            overlaps="-x"
        else
            overlaps=""
        fi
        echo "samtools mpileup -t DP,DP4 -L 1000 -d 1000 -m ${params.depth} \$anomolous \$BAQ \$overlaps -ugf ${params.ref} ${runname}/${name}.bam > ${runname}/tmp.mpileup" >> ${bashfile}
    fi

    echo "bcftools call -P ${params.prior} -O b -A -M -S ${runname}/${name}.ploidy -${params.call} ${runname}/tmp.mpileup > ${runname}/${name}.bcf" >> ${bashfile}
    echo "bcftools index ${runname}/${name}.bcf" >> ${bashfile}
    echo "bcftools call -P ${params.prior} -O b -A -M -v -S ${runname}/${name}.ploidy -${params.call} ${runname}/tmp.mpileup > ${runname}/${name}_variant.bcf" >> ${bashfile}
    echo "bcftools index ${runname}/${name}_variant.bcf" >> ${bashfile}

    if [ ${params.dirty} = "false" ]
    then
        echo "rm ${runname}/tmp.*" >> ${bashfile}
    fi

    if [ "${params.pseudosequence}" = "true" ]; then
        if [ "${params.call}" = "m" ]; then
            echo "bcf_2_pseudosequence.py -A -b ${runname}/${name}.bcf -B ${runname}/${name}.bam -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}" >> ${bashfile}
        elif [ "${params.call}" = "c" ]; then
            echo "bcf_2_pseudosequence.py -A -b ${runname}/${name}.bcf -B ${runname}/${name}.bam -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}" >> ${bashfile}
        fi
    fi
    """
}