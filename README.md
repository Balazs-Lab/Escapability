# HIV broadly neutralizing antibody escapability drives the therapeutic efficacy of vectored immunotherapy

This repository contains the following:
* Code used for analysis of viral escape from antibodies 
* Interactive viral escape data tables for each sample set, these can be used to explore the frequency of mutations at every site
* Code used for analysis of viral diversity prior to treatment with antibodies 


## System Requirements
The details about how to run the escapability pipeline can be found in the Viral Escape Analysis Directory in the Analysis_Notes.md file. The software required to run the escapability pipeline can be installed using conda and by manually installing CondonCaller software using python -m pip. While this software was developed and run on MacOS system, the pipline can be run on any standard machine that can run python and conda. 

The software required to run the escapability pipeline can be found in the environment file:

    viral_escape_environment.yml


## Installation Guide

To isntall the software required to run the escapability pipeline, ensure that the most recent version of Conda is installed on the computer. Additionally, the CodonCaller software will need to be downloaded from github and installed.  

Code for the codon calling software can be found in the following repository: https://github.com/Balazs-Lab/CodonCaller

The software required to run the escapability pipeline can be found in the environment file:

    viral_escape_environment.yml
     
The environment can be installed using conda

    conda env create -f viral_escape_environment.yml
    
Activate the conda environment and install the CodonCaller software 

    conda activate viral_escape
    cd /path/to/CodonCaller
    python -m pip install -e CodonCaller
     
    
## Demo / Instructions for use

A full demo to run the pipeline can be found in the Viral Escape Analysis Directory in the Analysis_Notes.md file, with the details of the pipeline found in the Snakefile. Once the conda environment and CondoCaller software are installed, the pipeline can be run using the snakemake command:

    snakemake --cores 4 

Specifying the number of cores will allow the processes to run in parallel. Total run time will depend on the number of samples, but a typical sample will finish the processing in less than an hour on a normal desktop computer. 

The output from the viral escape pipeline can be found in the Viral Escape Data directory.




