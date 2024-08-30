Viral Diversity Pipeline
========================

## Overview
This is the pipeline that was used to evaluate the diversity of HIV populations prior to anitbody treatment. Two toy raw data samples are included for demonstration purposes. If you are interested in re-running the entire analysis with all of the raw data, please contact us via the corresponding author's email address.

## System Requirements and Installation Guide
The software required to run the escapability pipeline can be installed using [Conda](https://conda.io). While this software was developed and run on a MacOS system, the pipeline can be run on any standard machine that can run Python and Conda. All other required software will be installed when creating the Conda environment:  

    conda env create -f environment.yml

## Demo

Before running the snakemake, be sure that all files have the correct naming convention. From the name, the computer should be able to determine the experiment number (i.e CD10), the mouse ID, the virus strain, and the week the sample was taken. Additionally, identical samples with multiple sequencing runs (multiple fastq files) can be merged into one file. 

The demo files are CD00-m000-00-JRCSF and CD00-m000-00-REJOc. Theses samples are two subsets of actual files which allow the user to test that the pipeline is integrated and running properly. To specify which samples are to be run in the pipeline, update the [Config.yaml](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Pipeline/config.yaml) file with the sample names and the fastq file path. 
  
Creat environment:

    conda env create -f environment.yml

Activate environment:

    conda activate ampseq
    
Make sure to build reference indexes:
  
  bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc-reference.fa
  bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF-reference.fa
  
  lofreq faidx data/reference_index/JRCSF-reference.fa
  lofreq faidx data/reference_index/REJOc-reference.fa


Run Snakemake:

    snakemake --cores 2
