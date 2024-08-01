# Multiple Mappings to BAM

## Clone the Repository

```
git clone git@github.com:LED-0102/multiple-mappings-to-bam.git
```

Move to dockerfiles directory

```
cd scripts/dockerfiles
```

Run the build script for convenience

```
chmod +x build_image.sh
./build_image.sh
```

Move back to root directory

```
cd ../../
```
Run the Nextflow script.

Replace the content between $ signs to the paths of the actual file
```
nextflow run scripts/main.nf --ref $Streptococcus_agalactiae_NGBS128_GCF_001552035_1.fa$ --program BWA --domapping True --human False --pairedend True --maxinsertsize 1000 --mininsertsize 50 --ssahaquality 30 --maprepeats False --GATK True --markdup True --detectOverlaps False --pseudosequence True --incref True --indels True --quality 50 --mapq 20 --depth 8 --stranddepth 3 --anomolous True --BAQ True --circular True --ratio 0.8 --prior 0.001 --call c --output Streptococcus_agalactiae_NGBS128_GCF_001552035_1_bwa --force False --filter 1 --tabfile False --alnfile False --raxml False --model GTRGAMMA --bootstrap 100 --keep False --LSF True --LSFQ normal --mem 5 --nodes 20 --dirty False --mapfiles $"31663_7#10_1.fastq.gz,31663_7#10_2.fastq.gz,31663_7#12_1.fastq.gz,31663_7#12_2.fastq.gz"$ -process.echo -resume
```

