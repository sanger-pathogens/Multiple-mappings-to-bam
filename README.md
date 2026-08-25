# multiple-mappings-to-bam

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[[_TOC_]]

## Pipeline overview

**multiple-mappings-to-bam** is a Nextflow DSL2 pipeline for mapping short-read paired-end sequencing data to a reference genome, calling variants, and optionally generating pseudosequences for phylogenetic analysis. It supports three mapping tools (BWA, SMALT, SSAHA2) and produces per-sample filtered BAM alignments and BCF variant calls.

The pipeline performs the following steps:

1. **Input** — paired FASTQ files are discovered from a read directory (see [Input](#input)).
2. **Mapping** — reads are aligned to the reference using BWA (default), SMALT, or SSAHA2.
3. **BAM processing** — the alignment is sorted; duplicates are optionally marked with Picard; GATK indel realignment is optionally applied; the BAM is filtered according to the selected filter mode.
4. **Variant calling** — Samtools mpileup generates per-position coverage and bcftools calls variants into a BCF.
5. **Pseudosequence generation** (optional, default: enabled) — a per-sample pseudosequence FASTA is derived from the BCF i.e. a consensus sequence is derived from the reference sequence where called variants or reference alleles are placed at their respective position and unknown bases `N` are used when no allele could be confidently called; indels are joined across samples; a multi-sample SNP summary is produced.

## Usage

### Quickstart

#### From source code

1. Clone this repository with its submodules:

   ```bash
   git clone --recurse-submodules https://github.com/sanger-pathogens/Multiple-mappings-to-bam.git
   cd multiple-mappings-to-bam
   ```

2. To run with `singularity`, use the `-profile singularity` option:

   ```bash
   nextflow run main.nf \
       -profile singularity \
       --read_dir /path/to/reads/ \
       --ref /path/to/reference.fasta \
       --outdir my_output
   ```

   :warning: If no profile is specified the pipeline will run with the Sanger HPC-specific configuration (see below).

3. Once the run has finished successfully and you have inspected the output, clean up intermediate files. The `work/` directory and `.nextflow.log` are useful for troubleshooting — do not delete them until you are satisfied the outputs are correct:

   ```bash
   rm -rf work .nextflow*
   ```

   Alternatively, use `nextflow clean` for more fine-grained control over which runs and intermediate files are removed.

#### Using on the Sanger 'farm' HPC

This pipeline is configured by default to run on the Sanger HPC under the `standard` profile (no need to specify it with `-profile`) so it can submit its tasks to the LSF scheduler to run as jobs on the HPC, and also to use singularity containers with adequate filesystem mountings and benefit from centralised singularity image caching.

First load the latest pipeline module:

```bash
module load multiple-mappings-to-bam
```

Then run on the command line with `multiple-mappings-to-bam <options>`. For instance, to see a help message:

```bash
multiple-mappings-to-bam --help
```

Submit to HPC job scheduler (LSF):

```bash
bsub -o output.o -e error.e -q oversubscribed -R "select[mem>4000] rusage[mem=4000]" -M4000 \
    multiple-mappings-to-bam \
        --read_dir /path/to/reads/ \
        --ref /path/to/reference.fasta \
        --outdir my_output
```

### Input

#### Read directory (`--read_dir`)

A path to a directory containing paired-end FASTQ files. Files must follow the naming convention `<sample_ID>_1.fastq.gz` / `<sample_ID>_2.fastq.gz`; the sample ID is derived from the filename prefix.

```
reads/
  sampleA_1.fastq.gz
  sampleA_2.fastq.gz
  sampleB_1.fastq.gz
  sampleB_2.fastq.gz
```

#### Reference (`--ref`)

Path to a reference FASTA file. Indexes are built automatically if not already present alongside the reference.

### Output

Results are written to `--outdir` (default: `./results`):

```
results/
  <sample_ID>_<program>/
    <sample_ID>.bcf              # BCF with all sites
    <sample_ID>_variant.bcf      # BCF with variant sites only
    <sample_ID>.ploidy           # Ploidy file used for variant calling
    <sample_ID>.mpileup          # Samtools mpileup output
    <sample_ID>_metrics.txt      # Picard duplicate-marking metrics (when --markdup true)
    <sample_ID>.mfa              # Per-sample pseudosequence FASTA (when --pseudosequence true)
  <ref>.aln                      # Multi-sample indel-joined pseudosequence alignment (when --pseudosequence true)
  <ref>.aln.out                  # Multi-sample SNP summary (when --pseudosequence true)
  <ref>.aln_summary.out          # SNP summary statistics (when --pseudosequence true)
```

### Parameters

**Input/output options**

| Option       | Type   | Default     | Description                                                              |
| ------------ | ------ | ----------- | ------------------------------------------------------------------------ |
| `--read_dir` | `path` | (required)  | Path to directory containing paired `*_1.fastq.gz`/`*_2.fastq.gz` files. |
| `--ref`      | `path` | (required)  | Path to the reference FASTA file.                                        |
| `--embl`     | `path` | `""`        | Path to reference annotation file (EMBL format).                         |
| `--outdir`   | `path` | `./results` | Directory where results are written.                                     |

---

**Mapping options**

| Option                 | Type      | Default | Description                                                                                                                                 |
| ---------------------- | --------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `--program`            | `string`  | `BWA`   | Mapping program. Options: `BWA`, `SMALT`, `SSAHA`.                                                                                          |
| `--human`              | `boolean` | `false` | Optimise SMALT k-mer size and step for mapping against the human genome (SMALT only).                                                       |
| `--maxinsertsize`      | `integer` | `1000`  | Maximum insert size for paired-end reads (SMALT/SSAHA only).                                                                                |
| `--mininsertsize`      | `integer` | `50`    | Minimum insert size for paired-end reads (SMALT/SSAHA only).                                                                                |
| `--ssahaquality`       | `integer` | `30`    | Minimum Phred base quality score (SSAHA only).                                                                                              |
| `--circular`           | `boolean` | `true`  | Treat contigs as circular (SSAHA only).                                                                                                     |
| `--allow_multimapping` | `boolean` | `false` | Map all reads including multi-mapping reads. Default excludes ambiguous mappings (SMALT only).                                              |
| `--nomapid`            | `float`   | `0`     | Minimum identity threshold for a mapping to be reported (SMALT only).                                                                       |
| `--GATK`               | `boolean` | `true`  | Run GATK indel realignment.                                                                                                                 |
| `--markdup`            | `boolean` | `true`  | Mark duplicate reads with Picard.                                                                                                           |
| `--detectOverlaps`     | `boolean` | `false` | Enable read-pair overlap detection.                                                                                                         |
| `--filter`             | `integer` | `1`     | BAM filtering mode: `1`=none, `2`=remove unmapped, `3`=properly paired only, `4`=split mapped/unmapped, `5`=split properly paired/unpaired. |
| `--publish_raw_bam`    | `boolean` | `false` | Publish raw bam files, as generated by mappers (and conversion to BAM format). Only filtered bam files will be output by default."          |

|

---

**Variant calling options**

| Option               | Type      | Default | Description                                                                                          |
| -------------------- | --------- | ------- | ---------------------------------------------------------------------------------------------------- |
| `--call`             | `string`  | `c`     | bcftools caller: `c`=consensus, `m`=multiallelic.                                                    |
| `--prior`            | `float`   | `0.001` | Prior probability that a site is non-reference. Higher values increase sensitivity to rare variants. |
| `--BAQ`              | `boolean` | `false` | Apply samtools base alignment quality (BAQ) recalibration.                                           |
| `--dontuseanomolous` | `boolean` | `false` | Exclude anomalous read pairs from mpileup.                                                           |

---

**Pseudosequence options**

| Option             | Type      | Default | Description                                                             |
| ------------------ | --------- | ------- | ----------------------------------------------------------------------- |
| `--pseudosequence` | `boolean` | `true`  | Generate pseudosequences from variant calls.                            |
| `--quality`        | `integer` | `50`    | Minimum base call quality for pseudosequence generation.                |
| `--mapq`           | `integer` | `20`    | Minimum mapping quality for pseudosequence generation.                  |
| `--depth`          | `integer` | `8`     | Minimum number of reads required to call a SNP.                         |
| `--stranddepth`    | `integer` | `3`     | Minimum number of reads per strand required to call a SNP.              |
| `--ratio`          | `float`   | `0.8`   | Minimum SNP/mapping quality ratio cutoff.                               |
| `--raxml`          | `boolean` | `false` | Run RAxML phylogeny on the pseudosequence alignment.                    |
| `--bootstrap`      | `integer` | `100`   | Number of RAxML bootstrap replicates. Set `0` to disable bootstrapping. |
| `--tabfile`        | `boolean` | `false` | Output a tab-delimited file of SNPs.                                    |
| `--alnfile`        | `boolean` | `false` | Output a SNP alignment file.                                            |

---

**Logging options**

| Option              | Type      | Default | Description                                            |
| ------------------- | --------- | ------- | ------------------------------------------------------ |
| `--monochrome_logs` | `boolean` | `false` | Output logs in plain ASCII (disable coloured logging). |

### Advanced usage

#### Choosing a mapping program

BWA is the default and recommended mapper for most use cases. Use SMALT for more control over repeat handling (`--allow_multimapping`) or for optimised human genome mapping (`--human`). Use SSAHA2 for circular contig support (`--circular`).

#### Disabling pseudosequence generation

To run mapping and variant calling only, without generating pseudosequences:

```bash
multiple-mappings-to-bam \
    --read_dir /path/to/reads/ \
    --ref /path/to/reference.fasta \
    --pseudosequence false \
    --outdir my_output
```

#### BAM filtering modes

The `--filter` option controls how the output BAM is filtered:

| Mode | Behaviour                                                   |
| ---- | ----------------------------------------------------------- |
| `1`  | No filtering (default)                                      |
| `2`  | Remove unmapped reads                                       |
| `3`  | Keep properly paired reads only                             |
| `4`  | Split mapped and unmapped reads into separate BAMs          |
| `5`  | Split properly paired and unpaired reads into separate BAMs |

## Dependencies

All dependencies are containerised. No external databases are required.

### Software versions

| Software | Version      | Image                                              |
| -------- | ------------ | -------------------------------------------------- |
| BWA      | 0.7.17-r1188 | `quay.io/ssd28/gsoc-experimental/run-bwa:0.0.2`    |
| SMALT    | 0.7.6        | `quay.io/ssd28/gsoc-experimental/run-smalt:0.0.2`  |
| SSAHA2   | 2.5.5        | `quay.io/sangerpathogens/ssaha2:v2.5.5_cv3`        |
| Samtools | 1.3          | `quay.io/ssd28/gsoc-experimental/samtools:1.3`     |
| Picard   | 1.126        | `quay.io/ssd28/gsoc-experimental/picard:1.126`     |
| GATK     | 3.7.0        | `quay.io/ssd28/gsoc-experimental/gatk:3.7.0`       |
| bcftools | 1.11         | `quay.io/ssd28/gsoc-experimental/bcftools:1.11-c1` |

See `modules/` for pinned container versions.

## Troubleshooting

- **No reads found**: ensure FASTQ files in `--read_dir` follow the `<sample_ID>_1.fastq.gz` / `<sample_ID>_2.fastq.gz` naming convention.
- **Reference index not found**: the pipeline builds BWA/SMALT/SSAHA indexes automatically. Ensure the reference directory is writable.
- **GATK indel realignment fails**: GATK 3.7 requires a sequence dictionary alongside the reference. This is generated automatically by the pipeline; ensure the reference directory is writable.
- **Resuming a failed run**: add `-resume` to your command to restart from cached intermediate results.
- For further help, check `.nextflow.log` and the per-process `.command.log` logs in the `work/` directory.

Sanger users may find [this page](https://ssg-confluence.internal.sanger.ac.uk/spaces/PaMI/pages/181078206/General+pipeline+info#Generalpipelineinfo-Troubleshootingafailedpipelinerunandsendingabugreport) useful for troubleshooting Nextflow pipeline runs.

## Credits

This pipeline was developed as part of the [Google Summer of Code 2024 program](https://summerofcode.withgoogle.com/archive/2024/projects/6g03n2ZD).

## Issues and Contributions

**GitHub users:** if you find an issue with this pipeline, or would like to suggest an improvement, please log an issue or open a pull request on this repository.

**Sanger users:** if you need internal support, you can raise an issue on the PAM Freshservice portal: https://sanger.freshservice.com/support/catalog/items/426
