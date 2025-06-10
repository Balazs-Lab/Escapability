# Viral Escape Pipeline

## Overview

This pipeline was used to identify HIV antibody escape mutants as profiled in the associated publication. The entire pipeline is included in this directory, along with a local copy of the [CodonCaller](https://github.com/adamn102/CodonCaller) Python package, which is **fully integrated** and does **not require separate installation**. All relevant Python scripts used in sample analysis are also included.

Two small raw data samples are provided for demonstration purposes. If you are interested in re-running the full analysis with the complete dataset, please contact the corresponding author. To view the final output of the pipeline directly, including compiled data tables, visit the [Viral Escape Data](https://github.com/Balazs-Lab/Escapability/tree/main/Viral%20Escape%20Data) repository.

## System Requirements and Installation Guide

The pipeline runs on any system that supports [Conda](https://conda.io) and Python. Although originally developed on macOS, it should work on any standard platform (Linux, Windows Subsystem for Linux, etc.).

To set up the environment:

```bash
conda env create -f environment.yml
conda activate viral_escape
```

**Note**: CodonCaller is now included locally in this repository and used directly by the pipeline. No manual installation is necessary.

## Demo

Before running the pipeline, ensure that all FASTQ files follow the correct naming convention. From the filename, the pipeline expects to infer:

* The experiment number (e.g., `CD10`)
* The mouse ID
* The virus strain
* The week the sample was taken

Identical samples with multiple sequencing runs (multiple FASTQ files) should be merged prior to analysis.

### Included Demo Data

The provided demo samples are:

* `CD00-m000-00-JRCSF`
* `CD00-m000-00-REJOc`

These are subsets of actual files that allow you to verify the pipeline is functioning properly.

To specify which samples to run, edit the [`config.yaml`](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/config.yaml) file with sample names and paths to FASTQ files.

### Running the Demo

1. **Create environment**:

   ```bash
   conda env create -f environment.yml
   ```

2. **Activate environment**:

   ```bash
   conda activate viral_escape
   ```

3. **Build reference indexes** (required once per reference):

   ```bash
   bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc-reference.fa
   bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF-reference.fa

   bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc.fa
   bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF.fa

   lofreq faidx data/reference_index/JRCSF-reference.fa
   lofreq faidx data/reference_index/REJOc-reference.fa
   ```

4. **Run the pipeline**:

   ```bash
   snakemake --cores 2
   ```
