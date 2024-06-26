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

process MAKEPILEUP_FROM_SAM {
    container 'quay.io/ssd28/gsoc-experimental/makepileup-from-sam:0.0.1'

    input:
    val runname
    val name
    val newsmalt
    val cmdline

    script:
    def version = getSamtoolsVersion()

    """

    #Detect version of SAMTOOLS, as "samtools sort" differs in usage depending on the version
    chevron=""
    suffix=""
    anomolous=""
    if echo "${version} > 1.2" | bc -l | grep -q 1; then
        chevron=">"
        suffix=".bam"
    fi

    #Sort and mark duplicates
    if [ ${params.markdup} ]
    then
        samtools sort ${runname}/tmp1.bam \$chevron ${runname}/tmpsort\$suffix
        picard MarkDuplicates INPUT=${runname}/tmpsort.bam OUTPUT=${runname}/tmp1.bam METRICS_FILE=${runname}/${name}_metrics.txt
        rm ${runname}/tmpsort.bam
    fi

    samtools sort ${runname}/tmp1.bam \$chevron ${runname}/${name}\$suffix
    samtools index ${runname}/${name}.bam
    rm ${runname}/tmp1.bam


    #Add read groups and fix smalt header
    #sed is used to substitute the SO:unknown to SO:coordinate in the header

    samtools view -H ${runname}/${name}.bam | sed 's/SO:unknown/SO:coordinate/g' | sed 's/\\\\x00//g' > ${runname}/tmphead.sam

    now=\$(date +'%Y-%m-%dT%H:%M:%S')

    if [ "${params.program}" = "smalt" ] || [ "${params.program}" = "SMALT" ]
    then
        echo \"@RG\\tID:${name}\\tCN:Sanger\\tDT:\$now\\tPG:SMALT\\tPL:ILLUMINA\\tSM:${name}\" >> ${runname}/tmphead.sam
        if [ ${params.domapping} ] && [ ${newsmalt} = "false" ]
        then
            smaltversion=\$( smalt version | grep Version | awk '{print \$2}' )

            #Unresolved issue below

            #echo "echo \"@PG\\tID:SMALT\\tPN:SMALT\\tCL:${cmdline}\\tVN:\$smaltversion\" >> ${runname}/tmphead.sam" >> ${bashfile}

            #Above line has unresolved issue

        
        fi
    elif [ "${params.program}" = "bwa" ] || [ "${params.program}" = "BWA" ]
    then
        echo '@RG\\tID:${name}\\tCN:Sanger\\tDT:\$now\\tPG:BWA MEM\\tPL:ILLUMINA\\tSM:${name}' >> ${runname}/tmphead.sam
    fi

    samtools view -b -o ${runname}/tmphead.bam -H ${runname}/${name}.bam
    samtools merge -c -p -f -r -h ${runname}/tmphead.sam ${runname}/tmp.bam ${runname}/${name}.bam ${runname}/tmphead.bam
    mv ${runname}/tmp.bam ${runname}/tmp1.bam
    rm ${runname}/${name}.bam

    #run GATK indel realignment if selected
    if [ ${params.GATK} ]
    then
        samtools index ${runname}/tmp1.bam
        cp ${params.ref} ${runname}/tmpref.fa
        samtools faidx ${runname}/tmpref.fa
        picard CreateSequenceDictionary R=${runname}/tmpref.fa O=${runname}/tmpref.dict
        gatk -I ${runname}/tmp1.bam -R ${runname}/tmpref.fa -T RealignerTargetCreator -o ${runname}/tmp.intervals
        gatk -I ${runname}/tmp1.bam -R ${runname}/tmpref.fa -T IndelRealigner --filter_bases_not_stored -targetIntervals ${runname}/tmp.intervals -o ${runname}/tmp.bam
        mv ${runname}/tmp.bam ${runname}/tmp1.bam
        rm ${runname}/tmp1.bam.bai ${runname}/tmpref.* ${runname}/tmp.intervals ${runname}/tmphead.*
    fi

    samtools sort ${runname}/tmp1.bam \$chevron ${runname}/tmp\$suffix
    rm ${runname}/tmp1.bam

    #Filter the bam file if requested
    if [ "${params.filter}" = "1" ]
    then
        mv ${runname}/tmp.bam ${runname}/${name}.bam
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

    #Index the bam file to get the bai file
    samtools index ${runname}/${name}.bam


    #produce the pileup file

    if [ ${params.anomolous} ]
    then
        anomolous=" -A "
    else
        anomolous=""
    fi


    #Make ploidy file for sample

    echo 'echo "${name}    1 > ${runname}/${name}.ploidy"' >> ${bashfile}

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
        samtools mpileup -t DP,DP4 -C 50 -L 1000 -d 1000 -m ${params.depth} \$anomolous \$BAQ \$overlaps -ugf ${params.ref} ${runname}/${name}.bam > ${runname}/tmp.mpileup
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
        samtools mpileup -t DP,DP4 -L 1000 -d 1000 -m ${params.depth} \$anomolous \$BAQ \$overlaps -ugf ${params.ref} ${runname}/${name}.bam > ${runname}/tmp.mpileup
    fi

    bcftools call -P ${params.prior} -O b -A -M -S ${runname}/${name}.ploidy -${params.call} ${runname}/tmp.mpileup > ${runname}/${name}.bcf
    bcftools index ${runname}/${name}.bcf
    bcftools call -P ${params.prior} -O b -A -M -v -S ${runname}/${name}.ploidy -${params.call} ${runname}/tmp.mpileup > ${runname}/${name}_variant.bcf
    bcftools index ${runname}/${name}_variant.bcf

    # clean up:

    if [ ${params.dirty} = "false" ]
    then
        rm ${runname}/tmp.*
    fi

    #Produce the pseudosequence if requested

    if [ "${params.pseudosequence}" = "true" ]; then
        if [ "${params.call}" = "m" ]; then
            bcf_2_pseudosequence.py -A -b ${runname}/${name}.bcf -B ${runname}/${name}.bam -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}
        elif [ "${params.call}" = "c" ]; then
            bcf_2_pseudosequence.py -A -b ${runname}/${name}.bcf -B ${runname}/${name}.bam -r ${params.ratio} -d ${params.depth} -D ${params.stranddepth} -q ${params.quality} -m ${params.mapq} -o ${runname}/${name}
        fi
    fi
    """
}
