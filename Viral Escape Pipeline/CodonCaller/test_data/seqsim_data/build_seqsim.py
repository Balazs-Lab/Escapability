from seqsim.seqsim import sequence

#direct to fileapth
test_data_fasta_filepath = 'reference_sequences/REJOc-reference.fa'

s = sequence(test_data_fasta_filepath)

#create fastq reference
s.generate_fastq_reference(301,301)

#load mutants
mutant_file_path = 'reference_sequences/REJOc-mutants.fa'
s.generate_mutant(mutant_file_path)

### WT Only Testing ###
s = sequence(test_data_fasta_filepath)
#create fastq reference
s.generate_fastq_reference(301,301)
#load mutants
mutant_file_path = 'reference_sequences/REJOc-mutants.fa'
s.generate_mutant(mutant_file_path)
#generate data
popsize = 10
WT = 1
s.generate_population(1,int(WT*popsize),s.fastq_reference)
output_file = 'output_files/WT_test'
s.write_fastq(output_file,read= "R1")
s.write_fastq(output_file,read= "R2")

### SNP Testing ###
s = sequence(test_data_fasta_filepath)
#create fastq reference
s.generate_fastq_reference(301,301)
#load mutants
mutant_file_path = 'reference_sequences/REJOc-mutants.fa'
s.generate_mutant(mutant_file_path)
#generate data
popsize = 10
REJOc_SNP = 1
WT = 1 - REJOc_SNP
s.generate_population(1,int(WT*popsize),s.fastq_reference)
s.generate_population(1,int(REJOc_SNP*popsize),s.mutants['REJOc_SNP'])
output_file = 'output_files/SNP_test'
s.write_fastq(output_file,read= "R1")
s.write_fastq(output_file,read= "R2")

### Insertion Testing ###
s = sequence(test_data_fasta_filepath)

#create fastq reference
s.generate_fastq_reference(301,301)
#load mutants
mutant_file_path = 'reference_sequences/REJOc-mutants.fa'
s.generate_mutant(mutant_file_path)
#generate data
popsize = 10
REJOc_IN = 1
WT = 1 - REJOc_IN
s.generate_population(1,int(WT*popsize),s.fastq_reference)
s.generate_population(1,int(REJOc_IN*popsize),s.mutants['REJOc_IN'])
output_file = 'output_files/IN_test'
s.write_fastq(output_file,read= "R1")
s.write_fastq(output_file,read= "R2")

### DEL Testing ###
s = sequence(test_data_fasta_filepath)
#create fastq reference
s.generate_fastq_reference(301,301)
#load mutants
mutant_file_path = 'reference_sequences/REJOc-mutants.fa'
s.generate_mutant(mutant_file_path)
#generate data
popsize = 10
REJOc_DEL = 1
WT = 1 - REJOc_DEL
s.generate_population(1,int(WT*popsize),s.fastq_reference)
s.generate_population(1,int(REJOc_DEL*popsize),s.mutants['REJOc_DEL'])
output_file = 'output_files/DEL_test'
s.write_fastq(output_file,read= "R1")
s.write_fastq(output_file,read= "R2")


### INDEL Testing ###
s = sequence(test_data_fasta_filepath)
#create fastq reference
s.generate_fastq_reference(301,301)
#load mutants
mutant_file_path = 'reference_sequences/REJOc-mutants.fa'
s.generate_mutant(mutant_file_path)
#generate data
popsize = 10
REJOc_INDEL = 1
WT = 1 - REJOc_INDEL
s.generate_population(1,int(WT*popsize),s.fastq_reference)
s.generate_population(1,int(REJOc_INDEL*popsize),s.mutants['REJOc_INDEL'])
output_file = 'output_files/INDEL_test'
s.write_fastq(output_file,read= "R1")
s.write_fastq(output_file,read= "R2")






## set population parmeters
#popsize = 10
#REJOc_IN = 0.1
#REJOc_DEL = 0.1
#REJOc_INDEL = 0.1
#REJOc_SNP = 0.1
#WT = 1 - (REJOc_IN + REJOc_DEL + REJOc_INDEL + REJOc_SNP)
#
#
#s.generate_population(1,int(WT*popsize),s.fastq_reference)
#s.generate_population(1,int(REJOc_IN*popsize),s.mutants['REJOc_IN'])
#s.generate_population(1,int(REJOc_DEL*popsize),s.mutants['REJOc_DEL'])
#s.generate_population(1,int(REJOc_INDEL*popsize),s.mutants['REJOc_INDEL'])
#s.generate_population(1,int(REJOc_SNP*popsize),s.mutants['REJOc_SNP'])
#
#output_file = 'output_files/testpop'
#s.write_fastq(output_file,read= "R1")
#s.write_fastq(output_file,read= "R2")
