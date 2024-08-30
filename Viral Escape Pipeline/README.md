Viral Escape Pipeline 
=====================

## Overview 
This is the pipeline that was used to identify the HIV antibody escape mutants profiled in the paper. The entire pipeline is included in this directory, including a local copy of the [CodonCaller](https://github.com/adamn102/CodonCaller) python package and all of the relevant python scripts used in sample analysis. Two toy raw data samples are included for demonstration purposes. If you are interested in re-running the entire analysis with all of the raw data, please contact us via the corresponding author's email address. If you would like to directly view the pipeline output data, please see the [Viral Escape Data](https://github.com/Balazs-Lab/Escapability/tree/main/Viral%20Escape%20Data) repository. Here the compiled data excel tables can be downloaded.

## System Requirements and Installation Guide
The software required to run the escapability pipeline can be installed using [Conda](https://conda.io) and manually installing CondonCaller software using Python -m pip. If you clone this directory, there is a copy of CodonCaller included, which can be used on the demo or provided data. While this software was developed and run on a MacOS system, the pipeline can be run on any standard machine that can run Python and Conda. All other required software will be installed when creating the Conda environment:  

    conda env create -f environment.yml

For the Escape Pipeline, once the initial environment is set up, the CodonCaller software will need to be installed. Activate the Conda environment and install the CodonCaller software:

    conda activate viral_escape
    cd /path/to/CodonCaller
    python -m pip install -e CodonCaller

The typical installation time for this software on a conventional computer is less than ten minutes.

## Demo

Before running the snakemake, be sure that all files have the correct naming convention. From the name, the computer should be able to determine the experiment number (i.e CD10), the mouse ID, the virus strain, and the week the sample was taken. Additionally, identical samples with multiple sequencing runs (multiple fastq files) can be merged into one file. 

The demo files are CD00-m000-00-JRCSF and CD00-m000-00-REJOc. Theses samples are two subsets of actual files which allow the user to test that the pipeline is integrated and running properly. To specify which samples are to be run in the pipeline, update the [Config.yaml](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/config.yaml) file with the sample names and the fastq file path.  

Create environment:

    conda env create -f environment.yml

Activate environment:

    conda activate viral_escape

Install CondonCaller:

     python -m pip install -e CodonCaller

Build reference indexes:

    bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc-reference.fa
    bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF-reference.fa
    
    bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc.fa
    bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF.fa
    
    lofreq faidx data/reference_index/JRCSF-reference.fa
    lofreq faidx data/reference_index/REJOc-reference.fa

Run Snakemake::

    snakemake --cores 2

