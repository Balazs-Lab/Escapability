Viral Escape Mutant Identification Data Pipeline  -
================================================

Analysis Notes:

Before running the snakemake, be sure that all files have the correct naming convention. From the name, the computer should be able to determine the experiment number (i.e CD10), the mouse ID, the virus strain, and the week the sample was taken. Additionally, identical samples with multiple sequencing runs (multiple fastq files) can be merged into one file. 

Additionally, make sure the computer environment has been created before processing samples. Use a Conda based environment, created from the environment.yml file located in the directory (we named the environment "viral_escape").

Finally, in order for the software to run, the CodonCaller package [https://github.com/adamn102/CodonCaller] needs to be installed. Be sure to do this within the viral_escape conda environment. 

Process Details:

create config file::

	ls data/raw_fastq_files/ | cut -d "_" -f 1 | uniq > names.txt

	while read p; do
	  ls data/raw_fastq_files/"$p"*R1*.gz | xargs cat > data/fastq_files/"$p"_R1.fastq.gz
	  ls data/raw_fastq_files/"$p"*R2*.gz | xargs cat > data/fastq_files/"$p"_R2.fastq.gz
	done <names.txt

	echo samples: > config.yaml
	awk '{ print "    " $0 ": data/fastq_files/" $0 "_R1.fastq.gz" }' < names.txt >> config.yaml

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

