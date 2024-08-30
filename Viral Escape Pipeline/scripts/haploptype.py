from CodonCaller.haplotype import *
import argparse
import pandas as pd


File = Haplotype(
    snakemake.input[2],
    snakemake.input[1],
    snakemake.input[0])
    
File.all_haplotypes
File.write_to_csv(snakemake.output[0])
    
