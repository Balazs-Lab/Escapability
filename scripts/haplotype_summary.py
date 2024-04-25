import pandas as pd
import csv
import numpy as np




def glycan_check(motif):
    glycan = True
    try:
        if motif[0] != "N":
            glycan = False
        if motif[1] == "P":
            glycan = False
        if motif[2] not in ["S","T"]:
            glycan = False
        if (len(motif) != 3):
            glycan = False
    except:
        glycan = False
    return glycan


### Analysis 1 - Summarize Haplotype Diversity by Sample

def average_haplotype(virus="REJOc",antibody="VRC07",min_freq = 0.00):
    sample_list = list(pd.read_csv("../data/metadata/"+ virus + "_"+ antibody + "_VL_figure_samples.csv")["id"])
    sample_list = [x for x in sample_list if str(x) != 'nan']

    df = pd.DataFrame()

        # read in sample haplotypes
    for sample in sample_list:
        
        experiment = sample.split("-")[0]
        mouse = sample.split("-")[1]
        week = sample.split("-")[2].split("k")[1:][0]
        
        sample_id  = experiment + "-" + mouse + "-" + week + "-" + virus
        haplotype_df = pd.read_csv("../data/haplotypes/" + experiment + "-" + mouse + "-" + week + "-" + virus + ".csv")
        haplotype_df['sample'] = sample_id
        #filter haplotype df
        haplotype_df = haplotype_df.loc[haplotype_df["count"] >= min_freq]
        df = pd.concat([df,haplotype_df])

    # compute the average frequency of every haplotype across the samples
    summary_df = df.fillna(0).groupby(["site","haplotype"])["count"].sum().reset_index().sort_values("site")

    summary_df['count'] = summary_df['count'] / len(set(df['sample']))
    summary_df[['site','reference']] = summary_df['site'].str.split('_', expand=True)
    summary_df['glycan'] = summary_df['haplotype'].map(glycan_check)
    summary_df['percent'] = summary_df['count']*100
    grouped_df = summary_df.groupby(['site','glycan'])['percent'].sum()

    #filter anything less than 1%
#    summary_df = summary_df.loc[summary_df["count"] >= min_freq]

    summary_df.to_csv("../data/escape_data/haplotypes/" + virus + "_" + antibody + "_haplotypes.csv")
    grouped_df.to_csv("../data/escape_data/haplotypes/" + virus + "_" + antibody + "_summary.csv")
        



#virus = 'JRCSF'
#antibody = 'control'
#min_freq = 0.00
#sample_list = list(pd.read_csv("../data/metadata/"+ virus + "_"+ antibody + "_VL_figure_samples.csv")["id"])
#sample_list = [x for x in sample_list if str(x) != 'nan']
#
#df = pd.DataFrame()
#
## read in sample haplotypes
#for sample in sample_list:
#    print(sample)
#    
#    experiment = sample.split("-")[0]
#    mouse = sample.split("-")[1]
#    week = sample.split("-")[2].split("k")[1:][0]
#    
#    sample_id  = experiment + "-" + mouse + "-" + week + "-" + virus
#    haplotype_df = pd.read_csv("../data/haplotypes/" + experiment + "-" + mouse + "-" + week + "-" + virus + ".csv")
#    haplotype_df['sample'] = sample_id
#    #filter haplotype df
#    haplotype_df = haplotype_df.loc[haplotype_df["count"] >= min_freq]
#    df = pd.concat([df,haplotype_df])
#
## compute the average frequency of every haplotype across the samples
#summary_df = df.groupby(["site","haplotype"])["count"].sum().reset_index().sort_values("site")
#
#summary_df['count'] = summary_df['count'] / len(set(df['sample']))
#summary_df[['site','reference']] = summary_df['site'].str.split('_', expand=True)
#summary_df['glycan'] = summary_df['haplotype'].map(glycan_check)
#grouped_df = summary_df.groupby(['site','glycan'])['count'].sum().reset_index()
#
#        
 

    
    


#
## Process data:
## Figure XX REJOc with Control
print("REJOc")
RC = average_haplotype(virus="REJOc",antibody="Control")
RV = average_haplotype(virus="REJOc",antibody="VRC07")
RP = average_haplotype(virus="REJOc",antibody="PGDM1400")
RN = average_haplotype(virus="REJOc",antibody="N6")

## Figure XX JRCSF with Control
print("JRCSF")
JC = average_haplotype(virus="JRCSF",antibody="Control")
JV = average_haplotype(virus="JRCSF",antibody="VRC07")
JP = average_haplotype(virus="JRCSF",antibody="PGDM1400")
JN = average_haplotype(virus="JRCSF",antibody="N6")

## Figure XX JRCSF with Control
print("JV")
JC = average_haplotype(virus="JV",antibody="Control")
JV = average_haplotype(virus="JV",antibody="VRC07")

print("JD")
JC = average_haplotype(virus="JD",antibody="Control")
JV = average_haplotype(virus="JD",antibody="VRC07")

print("JDV")
JC = average_haplotype(virus="JDV",antibody="Control")
JV = average_haplotype(virus="JDV",antibody="VRC07")

print("RV")
JC = average_haplotype(virus="RV",antibody="Control")
JV = average_haplotype(virus="RV",antibody="VRC07")

print("RD")
JC = average_haplotype(virus="RD",antibody="Control")
JV = average_haplotype(virus="RD",antibody="VRC07")

print("RDV")
JC = average_haplotype(virus="RDV",antibody="Control")
JV = average_haplotype(virus="RDV",antibody="VRC07")
