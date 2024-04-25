configfile: "config.yaml"

rule all:
    input:
        #expand("data/calls/{sample}.vcf", sample=config["samples"])
        #expand("data/var_freq/{sample}.csv", sample=config["samples"])
        #expand("data/sorted_reads/{sample}.bam.bai", sample=config["samples"])
        #expand("data/coverage/{sample}.csv", sample=config["samples"])
        expand("data/haplotypes/{sample}.csv", sample=config["samples"]),
        "data/var_freq_merged.csv",
        "data/coverage_merged.csv"
        


#filter sequences
rule filter_sequences:
    input:
        R1 = "data/fastq_files/{sample}_R1.fastq.gz",
        R2 = "data/fastq_files/{sample}_R2.fastq.gz"
    params:
        tm="0",
        qs="25", # the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
        ps="20" #  how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%
    log:
        "logs/fastq/{sample}.html"
    output:
        R1 = "data/fastq_filtered/{sample}_R1.fastq.gz",
        R2 = "data/fastq_filtered/{sample}_R2.fastq.gz"
    shell:
        "fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} "
        " --trim_tail1 {params.tm} --trim_tail2 {params.tm} --qualified_quality_phred {params.qs} "
        " --unqualified_percent_limit {params.ps} --html {log}"


def get_reference(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "data/reference_index/" + reference + "-reference.fa",

def get_bed(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "data/bed_file/" + reference + ".bed",
    
def get_primer_bed(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "data/bed_file/" + reference + "_primer.bed",
    
def get_gff(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "data/gff/" + reference + ".gff",

def get_refname(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "" + reference + "",
    
def get_cotig(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return reference,

def get_haplotype_bed(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "data/bed_file/" + reference + "_haplotypes.bed",

#bowtie2 mappping of filtered reads
rule bowtie2_map:
    input:
        #ref = "data/reference_index/Rp0-reference.fa" + reference,
        ref = get_reference,
        R1 = "data/fastq_filtered/{sample}_R1.fastq.gz",
        R2 = "data/fastq_filtered/{sample}_R2.fastq.gz"
    output:
        "data/mapped_reads/{sample}.bam"
    log:
        "logs/bowtie2/{sample}.log"
    threads: 8
    shell:
        "(bowtie2 -p {threads} -x {input.ref} -1 {input.R1} -2 {input.R2} "
        "-S {output}) 2> {log}"

rule samtools_sort:
    input:
        "data/mapped_reads/{sample}.bam"
    output:
        "data/sorted_reads/{sample}.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_index:
    input:
        "data/sorted_reads/{sample}.bam"
    output:
        "data/sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule call:
    input:
        fa=get_reference,
        bam="data/sorted_reads/{sample}.bam",
        bed=get_bed
    params:
        r=40,
        b=110,
        cds="Env",
    output:
        "data/calls/{sample}.csv",
        "data/coverage/{sample}.csv"
    script:
        "scripts/codon_caller.py"

rule merge_freq:
    input:
        expand("data/calls/{sample}.csv", sample=config["samples"])
    output:
        "data/var_freq_merged.csv"
    script:
        "scripts/merge_aa_freq.py"

rule merge_coverage:
    input:
        expand("data/coverage/{sample}.csv", sample=config["samples"])
    output:
        "data/coverage_merged.csv"
    script:
        "scripts/merge_coverage.py"

rule haplotype:
    input:
        fa=get_reference,
        bam="data/sorted_reads/{sample}.bam",
        bed=get_haplotype_bed
    output:
        "data/haplotypes/{sample}.csv",
    script:
        "scripts/haploptype.py"
