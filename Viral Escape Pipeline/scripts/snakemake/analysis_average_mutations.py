import pandas as pd
import csv
import numpy as np

### Overview ####
# This scrip will output summary data related to the average mutation frequency at each site
# and the types of amino acid changes at each site
# The output data will be in directories:
# escape_data/1-site_specific_mutation_summary
# escape_data/2-mutation_freq_summary


### Read input files
var_freq_df = pd.read_csv(snakemake.input["var_freq"])
coverage_df = pd.read_csv(snakemake.input["coverage"])
sample_list = list(pd.read_csv(snakemake.input["sample_list"])["id"])
exp_df = pd.read_csv(snakemake.input["exp_df"])


virus = snakemake.params["virus"]
antibody = snakemake.params["antibody"]
min_coverage = snakemake.params["min_coverage"]

alignment = pd.read_excel(snakemake.input["HXB2_alignment"],sheet_name=virus)



### Process Data Frames

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
# alignment = pd.read_excel("../../data/static/HXB2 Alignment.xlsx",sheet_name=virus)
position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)


summary_df["sum"] = summary_df[summary_df.columns[4:]].sum(axis=1)
summary_df["WT"] = 1-summary_df["sum"]

# write output file
outfile_path = snakemake.output["average_mutation"]
summary_df.to_csv(outfile_path)



