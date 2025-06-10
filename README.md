# HIV broadly neutralizing antibody escape dynamics drive the outcome of AAV vectored immunotherapy in humanized mice

## Pipelines and Data

This repository contains the following:

* [**Viral Escape Pipeline**](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/)
  Pipeline and code used for the analysis of viral escape from antibodies. This software was used to identify the primary antibody escape mutations.

* [**Viral Escape Data**](https://github.com/Balazs-Lab/Escapability/tree/main/Viral%20Escape%20Data)
  Output data from the escape pipeline in Excel format for each sample set. These tables allow exploration of the frequency of every amino acid at every site.

* [**Viral Diversity Pipeline**](https://github.com/Balazs-Lab/Escapability/tree/main/Viral%20Diversity%20Pipeline)
  Pipeline and code used for the analysis of viral diversity prior to antibody treatment.

* [**Escape Barrier Analysis**](https://github.com/Balazs-Lab/Escapability/tree/main/Escape%20Barrier%20Analysis)
  Data and scripts used for the escape barrier calculations.

## System Requirements

Each pipeline includes an `environment.yml` file listing software dependencies. This software was developed and tested on macOS but is fully compatible with any standard system that supports [Conda](https://conda.io) and Python (Linux, WSL, etc.).

* [Viral Escape Pipeline Requirements](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/environment.yml)
* [Viral Diversity Pipeline Requirements](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Pipeline/environment.yml)


## Demo / Instructions for Use

Full usage instructions are provided in the `README.md` file inside each pipeline directory:

* [Viral Escape Pipeline Demo](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/README.md)
* [Viral Diversity Pipeline Demo](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Pipeline/README.md)

Both pipelines include demo datasets and have been fully run on those examples, demonstrating all expected input/output files.

Each pipeline uses the [Snakemake](https://snakemake.readthedocs.io) workflow engine to manage execution. The core logic and rules are defined in:

* [Viral Escape Snakefile](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/Snakefile)
* [Viral Diversity Snakefile](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Diversity%20Pipeline/Snakefile)

Run time depends on sample number and compute resources.

## Viral Escape Pipeline Output Data

The processed [Viral Escape Data](https://github.com/Balazs-Lab/Escapability/tree/main/Viral%20Escape%20Data) used in the publication is available for download and further analysis.

After downloading the Excel file, navigate to the **Mutation Dashboard** tab to interact with the data. For each codon site in the HIV envelope, the following is available:

* Virus-strain-specific amino acid position alignment with HXB2
* Percent of total reads showing any amino acid change relative to wild type
* Breakdown of the types and frequencies of amino acid changes

## Escape Barrier Analysis

This analysis calculates the Escape Barrier Scores (AUC of the Escapability Plot), as described in the paper. In order to run this analysis, the full Viral Escape Pipeline will need to have been run and the Viral Haplotype Analysis will have need to have been performed.
