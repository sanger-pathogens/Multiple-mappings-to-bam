# Multiple Mappings to BAM

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[[_TOC_]]

## Introduction

This project is part of the Google Summer of Code 2024 program. It provides a Nextflow-based pipeline for processing DNA sequences and generating BAM files. The pipeline supports multiple mapping tools and various configurations to suit different research needs.

## Features

- Supports multiple mapping programs (BWA, SMALT, SSAHA)
- Handles paired-end and single-end reads
- Quality filtering and duplicate marking
- Generates pseudosequences
- Supports indel calling and variant detection
- Configurable parameters for advanced usage

### Getting started

### Running on the farm (Sanger HPC clusters)

1. Load nextflow and singularity modules:

   ```bash
   module load nextflow ISG/singularity
   ```

2. Either:

   - Clone this repository using `git clone --recurse-submodules`  
     OR
   - Use ready-made module: `module load Multiple-mappings-to-bam`  
     :warning: If using the ready-made module, replace `nextflow run main.nf` with `Multiple-mappings-to-bam` in all subsequent commands.

3. Start the pipeline  

   Example:


   ```bash
   nextflow run main.nf --manifest ./test_data/inputs/test_manifest.csv --outdir my_output
   ```

   It is good practice to submit a dedicated job for the nextflow master process (use the `oversubscribed` queue):

   ```bash
   bsub -o output.o -e error.e -q oversubscribed -R "select[mem>4000] rusage[mem=4000]" -M4000 nextflow run main.nf --ref absolute/path/to/ref/file --read_dir <path/to/directory/containing/reads> 
   ```

   See [usage](#usage) for all available pipeline options.

## Usage

```
Usage: 
    nextflow run main.nf 

Options:

-- outdir
    default: "./results"
    Output directory (optional)

-- read_dir 
    Absolute path to a directory containing reads (mandatory)

-- ref
    Absolute path to a directory containing reference DNA sequence (mandatory)

-- program
    default: bwa
    Mapping program you wish to include. Valid options [bwa, ssaha, smalt] (optional)

--domapping
    default: true
    Do not remap data (optional)

--human
    default: false
    Mapping against human (optional)
 
--pairedend
    default: true
    Reads are single ended (optional)

--maxinsertsize
    default: 1000
    Maximum insert size (optional)

--mininsertsize
    default: 50
    Minimum insert size (optional)

--ssahaquality
    default: 30
    Minimum ssaha quality score (optional)

--maprepeats
    default: false
    Randomly map repeats (optional)

--nomapid
    default: 0
    Minimum identity threshold (optional)

--GATK
    default: false
    Run GATK indel realignment (optional)

--markdup
    default: false
    Run Mark duplicates (optional)

--detectOverlaps
    default: false
    Enable read-pair overlap detection (optional)

--pseudosequence
    default: true
    Create pseudosequences (optional)

--incref
    default: true
    Include reference in pseudosequence alignment (optional)

--indels
    default: true
    Include small indels in pseudosequence alignment (optional)

--quality
    default: 50
    Minimum base call quality (optional)

--mapq
    default: 20
    Minimum mapping quality (optional)

--depth
    default: 8
    Minimum number of reads matching SNP (optional)

--stranddepth
    default: 3
    Minimum number of reads matching SNP per strand (optional)

--dontuseanomolous
    default: false
    Use anomolous reads in mpileup (optional)

--BAQ
    default: false
    Use samtools base alignment quality option (BAQ) (optional)

--circular
    default: true
    Contigs are circular (optional)

--ratio
    default: 0.8
    SNP/Mapping quality ratio cutoff (optional)

--prior
    default: 0.001
    Mutation rate (optional)

--call
    default: "c"
    bcftools caller (optional)

--embl
    default: ""
    Reference annotation (optional)

--output
    default: " "
    Output file prefix (optional)

--diroutput
    default: ""
    Output directory suffix (optional)

--force
    default: false
    Force overwrite of output files (optional)

--filter
    default: 1
    Filter or split bam file (optional)

--tabfile
    default: false
    Create tabfile of snps (optional)

--alnfile
    default: false
    Create snp alignment file (optional)

--raxml
    default: false
    Run phylogeny with RAxML (optional)

--model
    default: "GTRGAMMA"
    Model of evolution to use (optional)

--bootstrap
    default: 100
    Number of bootstrap replicates (optional)

--keep
    default: false
    If old mapping files are present, do not rerun them (optional)

--LSF
    default: true
    Use LSF to parallelise analyses (optional)

--LSFQ
    default: "normal"
    LSF queue to submit to (optional)

--mem
    default: 5
    Amount of memory required for analysis (Gb) (optional)

--nodes
    default: 20
    Maximum number of jobs to run on nodes in parallel (optional)

--dirty
    default: false
    Do not clean up temporary files (optional)

--docker
    default: false
    Run on docker containers (Default singularity) (optional)
```

## Support

For further information or help, don't hesitate to get in touch via [pam-informatics@sanger.ac.uk](mailto:pam-informatics@sanger.ac.uk).

