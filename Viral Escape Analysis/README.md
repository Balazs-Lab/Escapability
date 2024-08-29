Viral Escape Mutant Identification Data Pipeline  -
================================================

## System Requirements
The details about how to run the escapability pipeline can be found in the Viral Escape Analysis Directory in the Analysis_Notes.md file. The software required to run the escapability pipeline can be installed using conda and manually installing CondonCaller software using Python -m pip. If you clone this directory, there is a copy of CodonCaller included, which can be used on the demo or provided data. While this software was developed and run on a MacOS system, the pipeline can be run on any standard machine that can run Python and conda. 

## Demo

Before running the snakemake, be sure that all files have the correct naming convention. From the name, the computer should be able to determine the experiment number (i.e CD10), the mouse ID, the virus strain, and the week the sample was taken. Additionally, identical samples with multiple sequencing runs (multiple fastq files) can be merged into one file. 

The demo files are CD00-m000-00-JRCSF and CD00-m000-00-REJOc. Theses samples are two subsets of actual files which allow the user to test that the pipeline is itegrated and running properly. 

Additionally, make sure the computer environment has been created before processing samples. Use a Conda based environment, created from the environment.yml file located in the directory (we named the environment "viral_escape").

Finally, in order for the software to run, the [CodonCaller](https://github.com/adamn102/CodonCaller) package needs to be installed. Be sure to do this within the viral_escape conda environment. 

Creat environment::

    conda env create -f environment.yml

Activate environment::

    conda activate viral_escape

Install CondonCaller (if not already installed)::

     python -m pip install -e CodonCaller

Make sure to build reference indexes::

    bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc-reference.fa
    bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF-reference.fa
    
    bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc.fa
    bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF.fa
    
    lofreq faidx data/reference_index/JRCSF-reference.fa
    lofreq faidx data/reference_index/REJOc-reference.fa

Run Snakemake::

    snakemake --cores 2

