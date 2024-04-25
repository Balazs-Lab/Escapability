import pandas as pd
import csv
import numpy as np

### Analysis 1 - Codon Based Mapping of Average Mutations by Experimental Group

def average_mutaton(virus="REJOc",antibody="VRC07",min_coverage = 25):

    #read in variant frequency data
    var_freq_df = pd.read_csv("../data/var_freq_merged.csv")

    # read in coverage data
    coverage_df = pd.read_csv("../data/coverage_merged.csv")
    
    # read in sample list
    sample_list = list(pd.read_csv("../data/metadata/"+ virus + "_"+ antibody + "_VL_figure_samples.csv")["id"])

    # Read in metadata
    exp_df = pd.read_csv("../data/static/experiment_data.csv")

    var_freq_df['week'] = var_freq_df['week'].apply(lambda x: int(x.split('wk')[-1]))
    var_freq_df['id'] = var_freq_df['exp'] + "-" + var_freq_df['sample']
    var_freq_df = pd.merge(var_freq_df,exp_df[['id','antibody']],on='id',how='left')
    var_freq_df['id_antibody'] = var_freq_df['id'] + '-' + var_freq_df['antibody']

    coverage_df['week'] = coverage_df['week'].apply(lambda x: int(x.split('wk')[-1]))
    coverage_df['id'] = coverage_df['exp'] + "-" + coverage_df['sample']
    coverage_df = pd.merge(coverage_df,exp_df[['id','antibody']],on='id',how='left')
    coverage_df['id_antibody'] = coverage_df['id'] + '-' + coverage_df['antibody']


    escape_df = var_freq_df


    # list of final sequiencing time point for each mouse
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
#    print(aa_df)
    aa_df.sort_values("normalized_freq")


    aa_df.columns = ['aa_pos', 'aa_alt' ,'raw_freq','coverage','freq']

#    pivot_columns = ["aa_pos","*","A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    pivot_columns = ["aa_pos","*","A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","DEL","*i","Ai","Ri","Ni","Di","Ci","Ei","Qi","Gi","Hi","Ii","Li","Ki","Mi","Fi","Pi","Si","Ti","Wi","Yi","Vi"]
    pivot_df = aa_df.pivot(index=["aa_pos"],columns="aa_alt",values="freq").fillna(0).reset_index().sort_values("aa_pos").reindex(columns=pivot_columns).fillna(0)

    # align to HXB2 and ensure full range of samples
    alignment = pd.read_excel("../data/HXB2 Alignment.xlsx",sheet_name=virus)
    position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
    position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
    summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)

    summary_df["sum"] = summary_df[summary_df.columns[4:]].sum(axis=1)
    summary_df["WT"] = 1-summary_df["sum"]

    summary_df.to_csv("../data/escape_data/grouped_summary/" + virus + "_" + antibody + "_avg_freq.csv")
    
    out_df =summary_df[['HXB2_pos','sum']]
    
    out_df.columns = ["HXB2_pos",antibody]
    
    return out_df
            
   
#
## Process data:
## Figure XX REJOc with Control
RC = average_mutaton(virus="REJOc",antibody="Control")
RV = average_mutaton(virus="REJOc",antibody="VRC07")
RP = average_mutaton(virus="REJOc",antibody="PGDM1400")
RN = average_mutaton(virus="REJOc",antibody="N6")

# merge data
df = RC.merge(RV,how='left',on='HXB2_pos').merge(RP,how='left',on='HXB2_pos').merge(RN,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/REJOc Mutant Frequency Summary.csv")
#

## Figure 2C JRCSF with Control
JC = average_mutaton(virus="JRCSF",antibody="Control")
JV = average_mutaton(virus="JRCSF",antibody="VRC07")
JP = average_mutaton(virus="JRCSF",antibody="PGDM1400")
JN = average_mutaton(virus="JRCSF",antibody="N6")

# merge data
df = JC.merge(JV,how='left',on='HXB2_pos').merge(JP,how='left',on='HXB2_pos').merge(JN,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/JRCSF Mutant Frequency Summary.csv")

# Chimeras
control = average_mutaton(virus="R",antibody="Control")
antibody = average_mutaton(virus="R",antibody="VRC07")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/REJOc R Chimera Mutant Frequency Summary.csv")

control = average_mutaton(virus="RV",antibody="Control")
antibody = average_mutaton(virus="RV",antibody="VRC07")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/REJOc RV Chimera Mutant Frequency Summary.csv")

control = average_mutaton(virus="RD",antibody="Control")
antibody = average_mutaton(virus="RD",antibody="VRC07")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/REJOc RD Chimera Mutant Frequency Summary.csv")

control = average_mutaton(virus="RDV",antibody="Control")
antibody = average_mutaton(virus="RDV",antibody="VRC07")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/REJOc RDV Chimera Mutant Frequency Summary.csv")

control = average_mutaton(virus="J",antibody="Control")
antibody = average_mutaton(virus="J",antibody="VRC07")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/JRCSF J Chimera Mutant Frequency Summary.csv")

control = average_mutaton(virus="JV",antibody="Control")
antibody = average_mutaton(virus="JV",antibody="VRC07")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/JRCSF JV Chimera Mutant Frequency Summary.csv")

control = average_mutaton(virus="JD",antibody="Control")
antibody = average_mutaton(virus="JD",antibody="VRC07")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/JRCSF JD Chimera Mutant Frequency Summary.csv")

control = average_mutaton(virus="JDV",antibody="Control")
antibody = average_mutaton(virus="JDV",antibody="VRC07")
df = control.merge(antibody,how='left',on='HXB2_pos')
df.to_csv("../data/escape_data/JRCSF JDV Chimera Mutant Frequency Summary.csv")


### Analysis 2 - Site specific mutation frequency by sample

def mutation_freq_by_sample(virus="REJOc",antibody="VRC07",min_value=0.1,min_coverage=25):

    #read in variant frequency data
    var_freq_df = pd.read_csv("../data/var_freq_merged.csv")

    # read in coverage data
    coverage_df = pd.read_csv("../data/coverage_merged.csv")

    # read in sample list
    sample_list = list(pd.read_csv("../data/metadata/"+ virus + "_"+ antibody + "_VL_figure_samples.csv")["id"])

    # Read in metadata
    exp_df = pd.read_csv("../data/static/experiment_data.csv")

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

    sum_df = aa_df[['aa_pos','id','raw_freq']].groupby(by=["aa_pos","id"]).sum().reset_index()
    sum_df.columns = ["aa_pos","id","mut_freq"]
    pivot_df = sum_df.pivot(index=["aa_pos"],columns="id",values="mut_freq").fillna(0).reset_index().sort_values("aa_pos")

    # align to HXB2 and ensure full range of samples
    alignment = pd.read_excel("../data/HXB2 Alignment.xlsx",sheet_name=virus)
    position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
    position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
    summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)
    #summary_df = pivot_df.merge(position,how='right',on="aa_pos").fillna(0)


    summary_df.to_csv("../data/escape_data/sample_summary/" + virus + "_" + antibody + " " + str(min_value) +"_freq_by_sample.csv")

            
RC = mutation_freq_by_sample(virus="REJOc",antibody="Control",min_value=0)
RV = mutation_freq_by_sample(virus="REJOc",antibody="VRC07",min_value=0)
RP = mutation_freq_by_sample(virus="REJOc",antibody="PGDM1400",min_value=0)
RN = mutation_freq_by_sample(virus="REJOc",antibody="N6",min_value=0)
JC = mutation_freq_by_sample(virus="JRCSF",antibody="Control",min_value=0)
JV = mutation_freq_by_sample(virus="JRCSF",antibody="VRC07",min_value=0)
JP = mutation_freq_by_sample(virus="JRCSF",antibody="PGDM1400",min_value=0)
JN = mutation_freq_by_sample(virus="JRCSF",antibody="N6",min_value=0)

RC = mutation_freq_by_sample(virus="REJOc",antibody="Control",min_value=0.1)
RV = mutation_freq_by_sample(virus="REJOc",antibody="VRC07",min_value=0.1)
RP = mutation_freq_by_sample(virus="REJOc",antibody="PGDM1400",min_value=0.1)
RN = mutation_freq_by_sample(virus="REJOc",antibody="N6",min_value=0.1)
JC = mutation_freq_by_sample(virus="JRCSF",antibody="Control",min_value=0.1)
JV = mutation_freq_by_sample(virus="JRCSF",antibody="VRC07",min_value=0.1)
JP = mutation_freq_by_sample(virus="JRCSF",antibody="PGDM1400",min_value=0.1)
JN = mutation_freq_by_sample(virus="JRCSF",antibody="N6",min_value=0.1)

RV = mutation_freq_by_sample(virus="RV",antibody="VRC07",min_value=0.1)
RV = mutation_freq_by_sample(virus="RD",antibody="VRC07",min_value=0.1)
RDV = mutation_freq_by_sample(virus="RDV",antibody="VRC07",min_value=0.1)
RV = mutation_freq_by_sample(virus="RV",antibody="Control",min_value=0.1)
RV = mutation_freq_by_sample(virus="RD",antibody="Control",min_value=0.1)
RDV = mutation_freq_by_sample(virus="RDV",antibody="Control",min_value=0.1)

JV = mutation_freq_by_sample(virus="JV",antibody="VRC07",min_value=0.1)
JV = mutation_freq_by_sample(virus="JD",antibody="VRC07",min_value=0.1)
JDV = mutation_freq_by_sample(virus="JDV",antibody="VRC07",min_value=0.1)
JV = mutation_freq_by_sample(virus="JV",antibody="Control",min_value=0.1)
JV = mutation_freq_by_sample(virus="JD",antibody="Control",min_value=0.1)
JDV = mutation_freq_by_sample(virus="JDV",antibody="Control",min_value=0.1)


### Quality Check - Average Coverage per Sample
coverage_df = pd.read_csv("../data/coverage_merged.csv")
coverage_mean = coverage_df.groupby(['sample_id'])["COVERAGE"].mean().reset_index()
coverage_mean.columns = ["sample_id","mean_coverage"]
coverage_std = coverage_df.groupby(['sample_id'])["COVERAGE"].std().reset_index()
coverage_std.columns = ["sample_id","std_coverage"]
out_df = pd.merge(coverage_mean,coverage_std, how='outer',on='sample_id')
out_df.to_csv("../data/escape_data/coverage_summary.csv")



### Analysis 3 - Site specific mutations by sample

file_id = "CD5-m264-9"
virus = "REJOc"
def mutation_by_sample(file_id, virus):
    """
    Input: variant frequency file

    Output: HXB2 Aligned Site Specifc Mutation Frequencies (by AA change)
    also reports the coerage at each position
    """

    #read in coverage and variant data
    coverage_data = pd.read_csv("../data/coverage/"+file_id+"-"+virus+".csv")
    variant_data = pd.read_csv("../data/calls/"+file_id+"-"+virus+".csv")

    #merjge coverage data onto variant dataframe
    escape_df = variant_data[["POS_AA","REF_AA","ALT_AA","ALT_FREQ"]].merge(coverage_data[["POS_AA","COVERAGE"]],how='outer',on='POS_AA').fillna(0)
    escape_df.columns = ["aa_pos","REF_AA","ALT_AA","ALT_FREQ","COVERAGE"]

#    pivot_columns = ["aa_pos","COVERAGE","*","A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    pivot_columns = ["aa_pos","COVERAGE","*","A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","DEL","*i","Ai","Ri","Ni","Di","Ci","Ei","Qi","Gi","Hi","Ii","Li","Ki","Mi","Fi","Pi","Si","Ti","Wi","Yi","Vi"]

    pivot_df = escape_df.pivot(index=["aa_pos","COVERAGE"],columns="ALT_AA",values="ALT_FREQ").fillna(0).reset_index().sort_values("aa_pos").reindex(columns=pivot_columns).fillna(0)


    # align to HXB2 and ensure full range of samples
    alignment = pd.read_excel("../data/HXB2 Alignment.xlsx",sheet_name=virus)
    position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
    position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
    summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)
    #summary_df = pivot_df.merge(position,how='right',on="aa_pos").fillna(0)

    summary_df["sum"] = summary_df[summary_df.columns[5:]].sum(axis=1)
    summary_df["WT"] = 1-summary_df["sum"]


    summary_df.to_csv("../data/escape_data/sample_mutations/" + file_id + "_" + virus + "_mut_freq.csv")

# run on all the files var_freq
var_freq_df = pd.read_csv("../data/var_freq_merged.csv")
exp_df = pd.read_csv("../data/static/experiment_data.csv")

var_freq_df['week'] = var_freq_df['week'].apply(lambda x: int(x.split('wk')[-1]))
var_freq_df['id'] = var_freq_df['exp'] + "-" + var_freq_df['sample']
var_freq_df['file_id'] = var_freq_df['id']  + "-" + var_freq_df['week'].astype("string")
    
ids = list(exp_df['id'])
viruses = list(exp_df['virus'])
exp_dict = dict(zip(ids,viruses))

for file_id in list(set(var_freq_df['file_id'])):
    key = "-".join(file_id.split('-')[0:2])
    if key in exp_dict.keys():
        virus = exp_dict[key]
        mutation_by_sample(file_id,virus)

