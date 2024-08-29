Viral Diversity Data Pipeline  -
================================================

Analysis Notes:

Before running the snakemake, be sure that all files have the correct naming convention. From the name, the computer should be able to determine the experiment number (i.e CD10), the mouse ID, the virus strain, and the week the sample was taken. Additionally, identical samples with multiple sequencing runs (multiple fastq files) can be merged into one file. 

Additionally, make sure the computer environment have been created before processing samples.

Process Details:

  
Creat environment::

    conda env create -f environment.yml

Activate environment::

    conda activate ampseq
    
Make sure to build reference indexes:
  
  bowtie2-build data/reference_index/REJOc-reference.fa data/reference_index/REJOc-reference.fa
  bowtie2-build data/reference_index/JRCSF-reference.fa data/reference_index/JRCSF-reference.fa
  
  lofreq faidx data/reference_index/JRCSF-reference.fa
  lofreq faidx data/reference_index/REJOc-reference.fa


Run Snakemake::

    snakemake --cores 2
