# HIV broadly neutralizing antibody escapability drives the therapeutic efficacy of vectored immunotherapy

## Pipelines and Data

This repository contains the following:
* Viral Escape Analysis - pipeline and code used for the analysis of viral escape from antibodies, this was the software used to identify the primary antibody escape mutations 
* Viral Escape Data - viral escape Excel data tables for each sample set, these can be used to explore the frequency of mutations at every site
* Viral Diversity Analysis - pipeline and code used for analysis of viral diversity prior to treatment with antibodies 

## System Requirements
The details about how to run each pipeline can be found in their respective directories. 
Both pipelines contain an environment.yml file documenting the specific file dependencies. This software was developed and run on a MacOS system; however, the pipeline can be run on any standard machine that can run Python and [Conda](https://conda.io/).
 
[Viral Escape Pipeline Requirements](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Analysis/environment.yml)
[Viral Diversity Pipeline Requirements](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Analysis/environment.yml)  
 

## Installation Guide

To install the software required to run the escapability pipeline, ensure that the most recent version of [Conda](https://conda.io) is installed on the computer. Additionally, the [CodonCaller](https://github.com/Balazs-Lab/CodonCaller)  software will need to be downloaded from GitHub and installed. 

Once [Conda](https://conda.io) is installed, either pipeline environment can be downloaded using the following code:

    conda env create -f environment.yml

For the Escape Pipeline, once the initial environment is set up, the CodonCaller software will need to be installed. Activate the conda environment and install the CodonCaller software:

    conda activate viral_escape
    cd /path/to/CodonCaller
    python -m pip install -e CodonCaller

The typical installation time for this software on a conventional computer is less than ten minutes.
    

## Demo / Instructions for use

A full demo to run the pipeline can be found in the Viral Escape Analysis Directory in the Analysis_Notes.md file, with the details of the pipeline found in the Snakefile. Once the conda environment and CondoCaller software are installed, the pipeline can be run using the snakemake command:

    snakemake --cores 2 

Specifying the number of cores will allow the processes to run in parallel. Total run time will depend on the number of samples, but a typical sample will finish the processing in less than an hour on a normal desktop computer. 

The expected output from the viral escape pipeline is a .csv file calling the frequency of each mutation in a per mouse base that can be found in the Viral Escape Data directory. These files can later be merged together, generating a summary file with all data.

# Viral Escape Data

# Viral Diversity Analysis
