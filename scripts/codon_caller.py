from CodonCaller import *
import argparse

## parse arguments
#def main():
#    parser = argparse.ArgumentParser(description='Description of your tool')
#    parser.add_argument('bed', help='bed_file')
#    parser.add_argument('bam', help='bam_file')
#    parser.add_argument('cds', help='You MUST specify which CDS you want to call from')
#    parser.add_argument('ref', help='reference_fasta')
#    parser.add_argument('read_quality_threshold', help='read_quality_threshold')
#    parser.add_argument('base_quality_threshold', help='base_quality_threshold')
#    arser.add_argument('output', help='base_quality_threshold')
#
#
File = VariantCaller(
    snakemake.input[2],
    snakemake.input[1],
    snakemake.input[0],
    snakemake.params[0],
    snakemake.params[1])
    
File.process_sample(snakemake.params[2])
File.aa_count_to_freq(snakemake.params[2])
File.write_to_csv(snakemake.output[0])
File.write_coverage_to_csv(snakemake.output[1])

#vcf_to_ext_df(snakemake.input[0], snakemake.output[0])
