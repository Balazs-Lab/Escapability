Escape Barrier Analysis 
=======================

## Overview 

This analysis calculates the Escape Barrier Scores (AUC of the Escapability Plot), as described in the paper. In order to run this analysis, the full Viral Escape Pipeline will need to have been run and the Viral Haplotype Analysis will have need to have been performed. 

## Haplotype Analysis 
### Viral Escape Path haplotype categorization for each viral escape mouse sample.   
This directory contains the haplotype summary analysis of nonsynonymous mutations from the [Viral Escape Pipeline](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/). This data is passed into the Raw 0.1 tab and is used to determine haplotypes. Additional metadata is also incorporated, such as survival data, bNAb expression, post escape viral load (aka setpoint viral load). 
[Viral Escape Pipeline](https://github.com/Balazs-Lab/Escapability/blob/main/Viral%20Escape%20Pipeline/) must be run before this analysis can be performed.

There is one haplotype analysis sheet per bNAb - Virus combination presented in the paper. 


## Escape Barrier Score
### Viral Escape Path Escape Barrier Calculation and Data
This directory contains the escape barrier score Rscript, analysis input data, and analysis output data. 
### Input Data 
[Path Coordinates](https://github.com/Balazs-Lab/Escapability/blob/main/Escape%20Barrier%20Analysis/escape%20barrier%20score/data/Path%20Coordinates.xlsx): This table contains the compiled IC50 and Growth Rate (Doubling Time) data for each virus-antibody escape path combination. 
[Path Frequencies](https://github.com/Balazs-Lab/Escapability/blob/main/Escape%20Barrier%20Analysis/escape%20barrier%20score/data/Path%20Frequencies.xlsx): This table contains the haplotype specific frequency of each virus-antibody escape path combination. This data is used to scale the relative contributions of each esacpe path to the total Escape Barrier score.  

### Output Data


## Post Escape Viral Load 
