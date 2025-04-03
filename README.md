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
nextflow run main.nf --ref absolute/path/to/ref/file --program BWA --read_dir <path/to/directory/containing/reads>
````

### Directly running the pipeline

Without directly cloning the pipeline it can be run at the latest version like so

````sh
nextflow run sanger-pathogens/Multiple-mappings-to-bam --ref absolute/path/to/ref/file --program BWA --read_dir <path/to/directory/containing/reads>
````

The outputs would be generated in the directory specified by the `outdir` param (default `results/`)
