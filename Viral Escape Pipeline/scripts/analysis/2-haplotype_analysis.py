import pandas as pd
import csv
import numpy as np


### Analysis 2 - Site specific mutation frequency by sample - used to group samples for haplotype based analysis

def mutation_freq_by_sample(virus="REJOc",antibody="VRC07",min_value=0.1,min_coverage=25):

    # read in variant frequency data
    var_freq_df = pd.read_csv("../../data/7-merged_output/var_freq_merged.csv")

    # read in coverage data
    coverage_df = pd.read_csv("../../data/7-merged_output/coverage_merged.csv")

    # read in sample list
    sample_list = list(pd.read_csv("../../data/metadata/"+ virus + "_"+ antibody + "_VL_figure_samples.csv")["id"])

    # Read in metadata
    exp_df = pd.read_csv("../../data/static/experiment_data.csv")
    
    # organize metadata value
    var_freq_df['week'] = var_freq_df['week'].apply(lambda x: int(x.split('wk')[-1]))
    var_freq_df['id'] = var_freq_df['exp'] + "-" + var_freq_df['sample']
    var_freq_df = pd.merge(var_freq_df,exp_df[['id','antibody']],on='id',how='left')
    var_freq_df['id_antibody'] = var_freq_df['id'] + '-' + var_freq_df['antibody']

    coverage_df['week'] = coverage_df['week'].apply(lambda x: int(x.split('wk')[-1]))
    coverage_df['id'] = coverage_df['exp'] + "-" + coverage_df['sample']
    coverage_df = pd.merge(coverage_df,exp_df[['id','antibody']],on='id',how='left')
    coverage_df['id_antibody'] = coverage_df['id'] + '-' + coverage_df['antibody']

    escape_df = var_freq_df
    
    # filter by samples in sample list
    escape_df = escape_df[escape_df["sample_id"].isin(sample_list)].sort_values("POS_AA")
    coverage_df = coverage_df[coverage_df["sample_id"].isin(sample_list)]

    # filter my min coverage
    filtered_coverage = coverage_df[coverage_df['COVERAGE'] >= min_coverage]
    pass_filter = filtered_coverage.groupby(by=["POS_AA","id"]).size().reset_index()
    pass_filter.columns = ["POS_AA","id","sample_count"]

    # merge dfs and only keep positions that pass min coverage filter
    aa_df = escape_df.groupby(by=["POS_AA","ALT_AA","id"])["ALT_FREQ"].sum().reset_index().sort_values("POS_AA")
    aa_df = aa_df.merge(pass_filter,how='left',on=["POS_AA","id"]).fillna(0)

    # only keep values over the min frequency threshold
    aa_df = aa_df[aa_df.ALT_FREQ > min_value]

    aa_df.columns = ['aa_pos', 'aa_alt', 'id','raw_freq','coverage']

    pivot_columns = ["aa_pos","id","raw_freq","mut"]

    # create summary df for all mutations #
    
    sum_df = aa_df[['aa_pos','aa_alt','id','raw_freq']].groupby(by=["aa_pos","aa_alt","id"]).sum().reset_index()
    sum_df.columns = ["aa_pos","aa_alt","id","mut_freq"]
    pivot_df = sum_df.pivot(index=["aa_pos","aa_alt"],columns="id",values="mut_freq").fillna(0).reset_index().sort_values("aa_pos")

    # align to HXB2 and ensure full range of samples
    alignment = pd.read_excel("../../data/static/HXB2 Alignment.xlsx",sheet_name=virus)
    position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
    position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
    summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)
    #summary_df = pivot_df.merge(position,how='right',on="aa_pos").fillna(0)

    summary_df.to_csv("../../data/8-escape_data/3-haplotype_summary_all/" + virus + "_" + antibody + " " + str(min_value) +"_freq_by_sample.csv")
    
    # create summary df for only synonymous mutations #
    
    sum_df = aa_df[['aa_pos','id','raw_freq']].groupby(by=["aa_pos","id"]).sum().reset_index()
    sum_df.columns = ["aa_pos","id","mut_freq"]
    pivot_df = sum_df.pivot(index=["aa_pos"],columns="id",values="mut_freq").fillna(0).reset_index().sort_values("aa_pos")

    # align to HXB2 and ensure full range of samples
    alignment = pd.read_excel("../../data/static/HXB2 Alignment.xlsx",sheet_name=virus)
    position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
    position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
    summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)
    #summary_df = pivot_df.merge(position,how='right',on="aa_pos").fillna(0)

    summary_df.to_csv("../../data/8-escape_data/4-haplotype_summary_nonsyn/" + virus + "_" + antibody + " " + str(min_value) +"_freq_by_sample.csv")
    

print("REJO.c vs Control No Filter")
RC = mutation_freq_by_sample(virus="REJOc",antibody="Control",min_value=0)
print("REJO.c vs VRC07 No Filter")
RV = mutation_freq_by_sample(virus="REJOc",antibody="VRC07",min_value=0)
print("REJO.c vs PGDM1400 No Filter")
RP = mutation_freq_by_sample(virus="REJOc",antibody="PGDM1400",min_value=0)
print("REJO.c vs N6 No Filter")
RN = mutation_freq_by_sample(virus="REJOc",antibody="N6",min_value=0)

print("JR-CSF vs Control No Filter")
JC = mutation_freq_by_sample(virus="JRCSF",antibody="Control",min_value=0)
print("JR-CSF vs VRC07 No Filter")
JV = mutation_freq_by_sample(virus="JRCSF",antibody="VRC07",min_value=0)
print("JR-CSF vs PGDM1400 No Filter")
JP = mutation_freq_by_sample(virus="JRCSF",antibody="PGDM1400",min_value=0)
print("JR-CSF vs N6 No Filter")
JN = mutation_freq_by_sample(virus="JRCSF",antibody="N6",min_value=0)

print("REJO.c vs Control 10% Filter")
RC = mutation_freq_by_sample(virus="REJOc",antibody="Control",min_value=0.1)
print("REJO.c vs VRC07 10% Filter")
RV = mutation_freq_by_sample(virus="REJOc",antibody="VRC07",min_value=0.1)
print("REJO.c vs PGDM1400 10% Filter")
RP = mutation_freq_by_sample(virus="REJOc",antibody="PGDM1400",min_value=0.1)
print("REJO.c vs N6 10% Filter")
RN = mutation_freq_by_sample(virus="REJOc",antibody="N6",min_value=0.1)

print("JR-CSF vs Control 10% Filter")
JC = mutation_freq_by_sample(virus="JRCSF",antibody="Control",min_value=0.1)
print("JR-CSF vs VRC07 10% Filter")
JV = mutation_freq_by_sample(virus="JRCSF",antibody="VRC07",min_value=0.1)
print("JR-CSF vs PGDM1400 10% Filter")
JP = mutation_freq_by_sample(virus="JRCSF",antibody="PGDM1400",min_value=0.1)
print("JR-CSF vs N6 10% Filter")
JN = mutation_freq_by_sample(virus="JRCSF",antibody="N6",min_value=0.1)

print("REJO.c Chimera Experiment vs VRC07 10% Filter")
R = mutation_freq_by_sample(virus="R",antibody="VRC07",min_value=0.1)
print("RV Chimera vs VRC07 10% Filter")
RV = mutation_freq_by_sample(virus="RV",antibody="VRC07",min_value=0.1)
print("RD Chimera vs VRC07 10% Filter")
RD = mutation_freq_by_sample(virus="RD",antibody="VRC07",min_value=0.1)
print("RDV Chimera vs VRC07 10% Filter")
RDV = mutation_freq_by_sample(virus="RDV",antibody="VRC07",min_value=0.1)

print("REJO.c Chimera Experiment vs Control 10% Filter")
R = mutation_freq_by_sample(virus="R",antibody="Control",min_value=0.1)
print("RV Chimera vs Control 10% Filter")
RV = mutation_freq_by_sample(virus="RV",antibody="Control",min_value=0.1)
print("RD Chimera vs Control 10% Filter")
RD = mutation_freq_by_sample(virus="RD",antibody="Control",min_value=0.1)
print("RDV Chimera vs Control 10% Filter")
RDV = mutation_freq_by_sample(virus="RDV",antibody="Control",min_value=0.1)

print("JR-CSF Chimera Experiment vs VRC07 10% Filter")
J = mutation_freq_by_sample(virus="J",antibody="VRC07",min_value=0.1)
print("JV Chimera vs VRC07 10% Filter")
JV = mutation_freq_by_sample(virus="JV",antibody="VRC07",min_value=0.1)
print("JD Chimera vs VRC07 10% Filter")
JD = mutation_freq_by_sample(virus="JD",antibody="VRC07",min_value=0.1)
print("JDV Chimera vs VRC07 10% Filter")
JDV = mutation_freq_by_sample(virus="JDV",antibody="VRC07",min_value=0.1)

print("JR-CSF Chimera Experiment vs Control 10% Filter")
J = mutation_freq_by_sample(virus="J",antibody="Control",min_value=0.1)
print("JV Chimera vs Control 10% Filter")
JV = mutation_freq_by_sample(virus="JV",antibody="Control",min_value=0.1)
print("JD Chimera vs Control 10% Filter")
JD = mutation_freq_by_sample(virus="JD",antibody="Control",min_value=0.1)
print("JDV Chimera vs Control 10% Filter")
JDV = mutation_freq_by_sample(virus="JDV",antibody="Control",min_value=0.1)

