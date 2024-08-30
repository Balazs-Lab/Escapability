# HIV broadly neutralizing antibody escapability drives the therapeutic efficacy of vectored immunotherapy

## Pipelines and Data

This repository contains the following:
* [Viral Escape Pipeline](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/) - pipeline and code used for the analysis of viral escape from antibodies, this was the software used to identify the primary antibody escape mutations 
* [Viral Escape Data](https://github.com/Balazs-Lab/Escapability/tree/main/Viral%20Escape%20Data) - viral escape pipeline output data Excel tables for each sample set, these can be used to explore the frequency of amino acid at every site
* [Viral Diversity Pipeline](https://github.com/Balazs-Lab/Escapability/tree/main/Viral%20Diversity%20Pipeline) - pipeline and code used for analysis of viral diversity prior to treatment with antibodies 

## System Requirements
The details about how to run each pipeline can be found in their respective directories. 
Both pipelines contain an environment.yml file documenting the specific software dependencies. This software was developed and run on a MacOS system; however, the pipeline can be run on any standard machine that can run Python and [Conda](https://conda.io/).
 
[Viral Escape Pipeline Requirements](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/environment.yml)

[Viral Diversity Pipeline Requirements](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Pipeline/environment.yml)  
 

## Installation Guide

To install the software required to run the escapability pipeline, ensure that the most recent version of [Conda](https://conda.io) is installed on the computer. Additionally, to run the Viral Escape Pipeline the [CodonCaller](https://github.com/Balazs-Lab/CodonCaller) software will need to be downloaded from GitHub and installed. 

Once [Conda](https://conda.io) is installed, either pipeline environment can be downloaded using the following code:

    conda env create -f environment.yml

For the Escape Pipeline, once the initial environment is set up, the CodonCaller software will need to be installed. Activate the Conda environment and install the CodonCaller software:

    conda activate viral_escape
    cd /path/to/CodonCaller
    python -m pip install -e CodonCaller

The typical installation time for this software on a conventional computer is less than ten minutes.
    

## Demo / Instructions for use

A full demo on how to to run each pipeline can be found in the README.md file of each pipeline directory.

[Viral Escape Pipeline Demo](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/README.md)

[Viral Diversity Pipeline Demo](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Pipeline/README.md)

Currently, both pipelines are populated with example datasets that have been fully run on the pipeline, demonstrating the input and output files of the pipeline. 

Both pipelines use the [snakemake](https://snakemake.readthedocs.io) pipeline system to process the data, such that the structure and code of the pipeline live in each pipeline's respective snakemake file:

[Viral Escape Pipeline](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/Snakefile)

[Viral Diversity Pipeline](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Pipeline/Snakefile)  
 

Specifying the number of cores will allow the processes to run in parallel. Total run time will depend on the number of samples, but a typical sample will finish the processing in less than an hour on a normal desktop computer. 

The expected output from the viral escape pipeline is a .csv file calling the frequency of each mutation in a per mouse base that can be found in the Viral Escape Data directory. These files can later be merged together, generating a summary file with all data.

## Viral Escape Pipeline Output Data

The pooled and processed [Viral Escape Data](https://github.com/Balazs-Lab/Escapability/tree/main/Viral%20Escape%20Data) that was used in the paper can be downloaded for additional analysis. Once the Excel file is downloaded and opened, move to the Mutation Dashboard tab to directly interact with the data. Users can select any position across the HIV envelope to view the following on a per codon site basis:
* Virus strain specific amino acid position alignment with HXB2
* Percent of total reads with an amino acid change relative to wild type (the reference sequence)
* Breakdown of the various types of amino acid changes 
