#from CodonCaller import CodonCaller.variant_caller
import variant_caller 
import argparse


File = variant_caller.VariantCaller(
    snakemake.input[2],
    snakemake.input[1],
    snakemake.input[0],
    snakemake.params[0],
    snakemake.params[1])
    
File.process_sample(snakemake.params[2])
File.aa_count_to_freq(snakemake.params[2])
File.write_to_csv(snakemake.output[0])
File.write_coverage_to_csv(snakemake.output[1])

