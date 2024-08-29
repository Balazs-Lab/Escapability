# HIV broadly neutralizing antibody escapability drives the therapeutic efficacy of vectored immunotherapy

This repository contains the following:
* Code used for the analysis of viral escape from antibodies 
* Interactive viral escape data tables for each sample set, these can be used to explore the frequency of mutations at every site
* Code used for analysis of viral diversity prior to treatment with antibodies 


## System Requirements
The details about how to run the escapability pipeline can be found in the Viral Escape Analysis Directory in the Analysis_Notes.md file. The software required to run the escapability pipeline can be installed using conda and manually installing CondonCaller software using Python -m pip. If you clone this directory, there is a copy of CodonCaller included, which cna be used on the demo or provided data. While this software was developed and run on a MacOS system, the pipeline can be run on any standard machine that can run Python and conda. 

The software required to run the escapability pipeline can be found in the environment file:

    viral_escape_environment.yml


## Installation Guide

To install the software required to run the escapability pipeline, ensure that the most recent version of [Conda](https://conda.io) is installed on the computer. Additionally, the CodonCaller software will need to be downloaded from GitHub and installed.  

Code for the codon calling software can be found in the [CodonCaller](https://github.com/Balazs-Lab/CodonCaller) repository.

The software required to run the escapability pipeline can be found in the environment file:

    viral_escape_environment.yml
     
The environment can be installed using conda

    conda env create -f viral_escape_environment.yml
    
Activate the conda environment and install the CodonCaller software 

    conda activate viral_escape
    cd /path/to/CodonCaller
    python -m pip install -e CodonCaller

The typical installation time for this software on a conventional computer is less than ten minutes.
     
    
## Demo / Instructions for use

A full demo to run the pipeline can be found in the Viral Escape Analysis Directory in the Analysis_Notes.md file, with the details of the pipeline found in the Snakefile. Once the conda environment and CondoCaller software are installed, the pipeline can be run using the snakemake command:

    snakemake --cores 2 

Specifying the number of cores will allow the processes to run in parallel. Total run time will depend on the number of samples, but a typical sample will finish the processing in less than an hour on a normal desktop computer. 

The expected output from the viral escape pipeline is a .csv file calling the frequency of each mutation in a per mouse base that can be found in the Viral Escape Data directory. These files can later be merged together, generating a summary file with all data.
