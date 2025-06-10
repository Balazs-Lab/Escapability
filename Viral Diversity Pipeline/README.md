# Viral Diversity Pipeline

## Overview

This pipeline evaluates the genetic diversity of HIV populations prior to antibody treatment. Two toy raw data samples are included for demonstration purposes. If you're interested in running the full analysis on the complete dataset, please contact the corresponding author for access.

## System Requirements and Installation

The pipeline requires [Conda](https://conda.io) for environment and dependency management. While developed on macOS, it is compatible with any system that supports Python and Conda.

To set up the environment:

```bash
conda env create -f environment.yml
conda activate viral_diversity
```

All necessary tools (e.g., `snakemake`, `bowtie2`, `lofreq`, `samtools`, `fastp`) will be installed via the provided `environment.yml`.

## Demo Instructions

Before running the pipeline, ensure all filenames follow the expected naming convention. Each sample name should encode:

* Experiment ID (e.g., `CD00`)
* Mouse ID (e.g., `m000`)
* Virus strain (e.g., `JRCSF`)
* Collection time point (e.g., `00` for week 0)

If multiple sequencing runs exist for the same sample, merge the FASTQ files into one prior to execution.

### Included Demo Samples

Two demo samples are provided:

* `CD00-m000-00-JRCSF`
* `CD00-m000-00-REJOc`

These are subsets of actual samples used to verify that the pipeline is functioning correctly.

To specify which samples to process, update the `config.yaml` file with the sample names and corresponding FASTQ file paths:
[config.yaml](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Pipeline/config.yaml)

### Index Building (One-Time Setup)

Before running Snakemake, build the reference indices:

```bash
# Bowtie2 indices
bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc-reference.fa
bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF-reference.fa

# Lofreq faidx indexing
lofreq faidx data/reference_index/JRCSF-reference.fa
lofreq faidx data/reference_index/REJOc-reference.fa
```

### Running the Pipeline

Once the environment is activated and indices are built, run the full pipeline:

```bash
snakemake --cores 2
```
