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
- Optionally generates pseudosequences
- Supports indel calling and variant detection
- Configurable parameters for advanced usage

## Overview

### Mapping

Reads, paired or single-end, can be mapped to a reference by a selection of tools: BWA, SMALT or SSAHA. Choice of program enables various tool-specific options.

If mapping against the human genome with `--program SMALT`, using `--human` optimises kmer size and step size for faster, more memory-efficient mapping.

### Pileup

Optional BAM filtering has 5 modes to choose from (passed to the `--filter` option):

- `1` No filtering (Default)
- `2` Remove unmapped reads
- `3` Keep properly paired reads only
- `4` Split mapped and unmapped into separate bams
- `5` Split properly paired and unpaired reads into separate bams

### Pseudosequence Generation

Optional workflow that can be deactivated by setting `--pseudosequence false`.
RAxML phylogeny - can adjust bootstrapping (0 = none, 1-1000 to define bootstrap replicates)

## Getting started

### Running on the farm (Sanger HPC clusters)

1. Load nextflow and singularity modules:

   ```bash
   module load nextflow ISG/singularity
   ```

2. Either:

   - Clone this repository using `git clone --recurse-submodules`  
     OR
   - Use ready-made module: `module load multiple-mappings-to-bam`  
     :warning: If using the ready-made module, replace `nextflow run main.nf` with `multiple-mappings-to-bam` in all subsequent commands.

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
    nextflow run main.nf [options]

Options:

Input/output options

-- read_dir
    Absolute path to a directory containing reads (mandatory)

-- ref
    Absolute path to a directory containing reference DNA sequence (mandatory)

-- embl
    default: ""
    Reference annotation

-- outdir
    default: "./results"
    Output directory

-- diroutput
    default: ""
    Output directory suffix

-- output
    default: " "
    Output file prefix

Mapping options


-- program
    default: BWA
    Mapping program you wish to include: BWA|SSAHA|SMALT

-- domapping
    default: true
    Do not remap data

-- human
    default: false
    Optimise SMALT for mapping against human genome

-- pairedend
    default: true
    Set to false for single-end reads

-- maxinsertsize
    default: 1000
    Maximum insert size for paired-end reads (SMALT/SSAHA only)

-- mininsertsize
    default: 50
    Minimum insert size for paired-end reads (SMALT/SSAHA only)

-- ssahaquality
    default: 30
    Minimum Phred base quality score (SSAHA only)

-- circular
    default: true
    Contigs are circular (SSAHA only)

-- maprepeats
    default: false
    Map all reads, including repeats (even ambiguous mappings). Default is false to exclude multi-mapping reads. (SMALT only)

-- nomapid
    default: 0
    Minimum identity threshold, as a float, for mapping to be reported (SMALT only)

-- GATK
    default: false
    Run GATK indel realignment (optional)

-- markdup
    default: true
    Mark duplicates with Picard (optional)

-- detectOverlaps
    default: false
    Enable read-pair overlap detection (optional)

-- filter
    default: 1
    Filtering mode for bam file (1=No filter, 2=remove unmapped, 3=properly paired only, 4=split mapped/unmapped, 5=split properly paired/unpaired)

Variant calling options:

-- call
    default: c
    bcftools caller (c=consensus, m=multiallelic)

-- prior
    default: 0.001
    Sets the prior probability that a site is non-reference for variant calling, higher values increase sensitivity to rare variants (optional)

-- BAQ
    default: false
    Use samtools base alignment quality option (BAQ) (optional)

-- dontuseanomolous
    default: false
    Use anomalous reads in mpileup (optional)

Pseudosequence options:

-- pseudosequence
    default: true
    Create pseudosequences (optional)

-- incref
    default: true
    Include reference in pseudosequence alignment (optional)

-- indels
    default: true
    Include small indels in pseudosequence alignment (optional)

-- quality
    default: 50
    Minimum base call quality (optional)

-- mapq
    default: 20
    Minimum mapping quality (optional)

-- depth
    default: 8
    Minimum number of reads matching SNP (optional)

-- stranddepth
    default: 3
    Minimum number of reads matching SNP per strand (optional)

-- ratio
    default: 0.8
    SNP/Mapping quality ratio cutoff (optional)

-- raxml
    default: false
    Run phylogeny with RAxML (optional)

-- model
    default: GTRGAMMA
    Model for RAxML: GTRCAT|GTRMIX|GTRGAMMA (optional)

-- bootstrap
    default: 100
    Number of bootstrap replicates for RAxML, 0=No bootstrap (optional)

-- tabfile
    default: false
    Create tabfile of snps (optional)

-- alnfile
    default: false
    Create snp alignment file (optional)

Job submission and workflow options

-- LSF
    default: true
    Use LSF to parallelise analyses (optional)

-- LSFQ
    default: "normal"
    LSF queue to submit to (optional)

-- mem
    default: 5
    Amount of memory required for analysis (Gb) (optional)

-- nodes
    default: 20
    Maximum number of jobs to run on nodes in parallel (optional)

-- force
    default: false
    Force overwrite of output files (optional)

-- dirty
    default: false
    Do not clean up temporary files (optional)
```

## Support

For further information or help, don't hesitate to get in touch via [pam-informatics@sanger.ac.uk](mailto:pam-informatics@sanger.ac.uk).
