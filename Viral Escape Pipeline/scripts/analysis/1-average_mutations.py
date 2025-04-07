import pandas as pd
import csv
import numpy as np

### Overview ####
# This scrip will output summary data related to the average mutation frequency at each site
# and the types of amino acid changes at each site
# The output data will be in directories:
# escape_data/1-site_specific_mutation_summary
# escape_data/2-mutation_freq_summary


### Analysis 1 - Codon Based Mapping of Average Mutations by Experimental Group

def average_mutaton(virus="REJOc",antibody="VRC07",min_coverage = 25):

    #read in variant frequency data
    var_freq_df = pd.read_csv("../../data/7-merged_output/var_freq_merged.csv")

    # read in coverage data
    coverage_df = pd.read_csv("../../data/7-merged_output/coverage_merged.csv")
    
    # read in sample list
    sample_list = list(pd.read_csv("../../data/metadata/"+ virus + "_"+ antibody + "_VL_figure_samples.csv")["id"])

    # Read in metadata
    exp_df = pd.read_csv("../../data/static/experiment_data.csv")

    var_freq_df['week'] = var_freq_df['week'].apply(lambda x: int(x.split('wk')[-1]))
    var_freq_df['id'] = var_freq_df['exp'] + "-" + var_freq_df['sample']
    var_freq_df = pd.merge(var_freq_df,exp_df[['id','antibody']],on='id',how='left')
    var_freq_df['id_antibody'] = var_freq_df['id'] + '-' + var_freq_df['antibody']

    coverage_df['week'] = coverage_df['week'].apply(lambda x: int(x.split('wk')[-1]))
    coverage_df['id'] = coverage_df['exp'] + "-" + coverage_df['sample']
    coverage_df = pd.merge(coverage_df,exp_df[['id','antibody']],on='id',how='left')
    coverage_df['id_antibody'] = coverage_df['id'] + '-' + coverage_df['antibody']


    escape_df = var_freq_df


    # list of final sequencing time point for each mouse
    meta_df = escape_df[["exp","sample","week","sample_id"]].drop_duplicates().reset_index()
    max_week_indices = meta_df.groupby("sample")["week"].idxmax()
    max_week_df = meta_df.loc[max_week_indices].reset_index()[["exp","sample","week","sample_id"]]

    # filter by samples in sample list
    escape_df = escape_df[escape_df["sample_id"].isin(sample_list)].sort_values("POS_AA")
    coverage_df = coverage_df[coverage_df["sample_id"].isin(sample_list)]

    # count unique samples in each group by coverage
    filtered_coverage = coverage_df[coverage_df['COVERAGE'] >= min_coverage]

    pass_filter = filtered_coverage.groupby(by=["POS_AA"]).size().reset_index()
    pass_filter.columns = ["POS_AA","sample_count"]
    
    # merge dfs and normalize by coverage
    aa_df = escape_df.groupby(by=["POS_AA","ALT_AA"])["ALT_FREQ"].sum().reset_index().sort_values("POS_AA")
    aa_df = aa_df.merge(pass_filter,how='left',on="POS_AA").fillna(0)
    aa_df["normalized_freq"] = aa_df['ALT_FREQ'] / aa_df["sample_count"]

    aa_df.sort_values("normalized_freq")

    aa_df.columns = ['aa_pos', 'aa_alt' ,'raw_freq','coverage','freq']

    pivot_columns = ["aa_pos","*","A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","DEL","*i","Ai","Ri","Ni","Di","Ci","Ei","Qi","Gi","Hi","Ii","Li","Ki","Mi","Fi","Pi","Si","Ti","Wi","Yi","Vi"]
    pivot_df = aa_df.pivot(index=["aa_pos"],columns="aa_alt",values="freq").fillna(0).reset_index().sort_values("aa_pos").reindex(columns=pivot_columns).fillna(0)

    # align to HXB2 and ensure full range of samples
    alignment = pd.read_excel("../../data/static/HXB2 Alignment.xlsx",sheet_name=virus)
    position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
    position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
    summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)


    summary_df["sum"] = summary_df[summary_df.columns[4:]].sum(axis=1)
    summary_df["WT"] = 1-summary_df["sum"]

    summary_df.to_csv("../../data/8-escape_data/2-mutation_freq_summary/" + virus + "_" + antibody + "_avg_freq.csv")
    
    out_df =summary_df[['HXB2_pos','sum']]
    
    out_df.columns = ["HXB2_pos",antibody]
    
    return out_df
            

## Process data:
## REJOc with Control
print("REJO.c vs Control")
RC = average_mutaton(virus="REJOc",antibody="Control")

print("REJO.c vs VRC07")
RV = average_mutaton(virus="REJOc",antibody="VRC07")

print("REJO.c vs PGDM1400")
RP = average_mutaton(virus="REJOc",antibody="PGDM1400")

print("REJO.c vs N6")
RN = average_mutaton(virus="REJOc",antibody="N6")

# merge and write data
print("REJO.c Summary")
df = RC.merge(RV,how='left',on='HXB2_pos').merge(RP,how='left',on='HXB2_pos').merge(RN,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/REJOc Mutant Frequency Summary.csv")


## JRCSF with Control
print("JR-CSF vs Control")
JC = average_mutaton(virus="JRCSF",antibody="Control")

print("JR-CSF vs VRC07")
JV = average_mutaton(virus="JRCSF",antibody="VRC07")

print("JR-CSF vs PGDM1400")
JP = average_mutaton(virus="JRCSF",antibody="PGDM1400")

print("JR-CSF vs N6")
JN = average_mutaton(virus="JRCSF",antibody="N6")

# merge and write data
print("JR-CSF Summary")
df = JC.merge(JV,how='left',on='HXB2_pos').merge(JP,how='left',on='HXB2_pos').merge(JN,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/JRCSF Mutant Frequency Summary.csv")

## Chimeras
# REJO.c
print("REJO.c Chimera Experiment vs Control")
control = average_mutaton(virus="R",antibody="Control")

print("REJO.c Chimera Experiment vs VRC07")
antibody = average_mutaton(virus="R",antibody="VRC07")

print("REJO.c Chimera Experiment Summary")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/REJOc R Chimera Mutant Frequency Summary.csv")

# RV
print("RV Chimera vs Control")
control = average_mutaton(virus="RV",antibody="Control")

print("RV Chimera vs VRC07")
antibody = average_mutaton(virus="RV",antibody="VRC07")

print("RV Chimera Summary")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/REJOc RV Chimera Mutant Frequency Summary.csv")

# RD
print("RD Chimera vs Control")
control = average_mutaton(virus="RD",antibody="Control")

print("RD Chimera vs VRC07")
antibody = average_mutaton(virus="RD",antibody="VRC07")

print("RD Chimera Summary")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/REJOc RD Chimera Mutant Frequency Summary.csv")

# RDV
print("RDV Chimera vs Control")
control = average_mutaton(virus="RDV",antibody="Control")

print("RDV Chimera vs VRC07")
antibody = average_mutaton(virus="RDV",antibody="VRC07")

print("RDV Chimera Summary")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/REJOc RDV Chimera Mutant Frequency Summary.csv")

# JR-CSF
print("JR-CSF Chimera Experiment vs Control")
control = average_mutaton(virus="J",antibody="Control")

print("JR-CSF Chimera Experiment vs VRC07")
antibody = average_mutaton(virus="J",antibody="VRC07")

print("JR-CSF Chimera Experiment Summary")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/JRCSF J Chimera Mutant Frequency Summary.csv")

# JV
print("JV Chimera vs Control")
control = average_mutaton(virus="JV",antibody="Control")

print("JV Chimera vs VRC07")
antibody = average_mutaton(virus="JV",antibody="VRC07")

print("JV Chimera Summary")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/JRCSF JV Chimera Mutant Frequency Summary.csv")

# JD
print("JD Chimera vs Control")
control = average_mutaton(virus="JD",antibody="Control")

print("JD Chimera vs VRC07")
antibody = average_mutaton(virus="JD",antibody="VRC07")

print("JD Chimera Summary")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/JRCSF JD Chimera Mutant Frequency Summary.csv")

# JDV
print("JDV Chimera vs Control")
control = average_mutaton(virus="JDV",antibody="Control")

print("JDV Chimera vs VRC07")
antibody = average_mutaton(virus="JDV",antibody="VRC07")

print("JDV Sumnmary")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../../data/8-escape_data/1-site_specific_mutation_summary/JRCSF JDV Chimera Mutant Frequency Summary.csv")
