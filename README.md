# Multiple Mappings to BAM

## Description

This project is part of the Google Summer of Code 2024 program. It provides a Nextflow-based pipeline for processing DNA sequences and generating BAM files. The pipeline supports multiple mapping tools and various configurations to suit different research needs.

## Features

- Supports multiple mapping programs (e.g., BWA, SMALT, SSAHA)
- Handles paired-end and single-end reads
- Quality filtering and duplicate marking
- Generates pseudosequences
- Supports indel calling and variant detection
- Configurable parameters for advanced usage

## Installation

To install the necessary dependencies, ensure you have Nextflow installed. You can refer to the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html) for detailed instructions.


````sh
curl -s https://get.nextflow.io | bash
````
Additionally, you need to install Docker or Singularity. Follow the respective installation guides:

[Docker Installation Guide](https://docs.docker.com/engine/install/)  
[Singularity Installation Guide](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)  

## Usage

Clone the repository and navigate to the project directory:

````sh
git clone
cd multiple-mappings-to-bam
````

### Using the input CLI

Use the input CLI to pass the parameters. Run the CLI:

````sh
python input/input.py
````

For detailed instructions on how to use the input CLI, refer to the USAGE.md file.

### Directly running the pipeline

Instead of using the input CLI, you can run the pipeline directly using the nextflow run command. It is not recommended as some input validation is handled by the CLI.  
Here's an example usage:

````sh
nextflow run scripts/main.nf --ref absolute/path/to/ref/file --program BWA --domapping True --human False --pairedend True --maxinsertsize 1000 --mininsertsize 50 --ssahaquality 30 --maprepeats False --GATK True --markdup True --detectOverlaps False --pseudosequence True --incref True --indels True --quality 50 --mapq 20 --depth 8 --stranddepth 3 --anomolous True --BAQ True --circular True --ratio 0.8 --prior 0.001 --call c --output Streptococcus_agalactiae_NGBS128_GCF_001552035_1_bwa --force False --filter 1 --tabfile False --alnfile False --raxml False --model GTRGAMMA --bootstrap 100 --keep False --LSF True --LSFQ normal --mem 5 --nodes 20 --dirty False --mapfiles "absolute/path/to/ref/file1,absolute/path/to/ref/file2,..(Add more files)" -process.echo -resume
````

The outputs would be generated in the directory specified by the `outdir` param (default `results/`)