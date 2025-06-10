import pandas as pd
import csv
import numpy as np



### Read input files
var_freq_df = pd.read_csv(snakemake.input["var_freq"])
coverage_df = pd.read_csv(snakemake.input["coverage"])
sample_list = list(pd.read_csv(snakemake.input["sample_list"])["id"])
exp_df = pd.read_csv(snakemake.input["exp_df"])


virus = snakemake.params["virus"]
antibody = snakemake.params["antibody"]
min_coverage = snakemake.params["min_coverage"]
min_value = snakemake.params["min_value"]

alignment = pd.read_excel(snakemake.input["HXB2_alignment"],sheet_name=virus)


### Analysis 2 - Site specific mutation frequency by sample - used to group samples for haplotype based analysis

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

# write output file
outfile_path = snakemake.output["haplotypes"]
summary_df.to_csv(outfile_path)

summary_df.to_csv("../../data/8-escape_data/3-haplotype_summary_all/" + virus + "_" + antibody + " " + str(min_value) +"_freq_by_sample.csv")

# create summary df for only synonymous mutations #

sum_df = aa_df[['aa_pos','id','raw_freq']].groupby(by=["aa_pos","id"]).sum().reset_index()
sum_df.columns = ["aa_pos","id","mut_freq"]
pivot_df = sum_df.pivot(index=["aa_pos"],columns="id",values="mut_freq").fillna(0).reset_index().sort_values("aa_pos")

# align to HXB2 and ensure full range of samples
position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)
#summary_df = pivot_df.merge(position,how='right',on="aa_pos").fillna(0)

# write output file
outfile_path = snakemake.output["haplotypes"]
summary_df.to_csv(outfile_path)

summary_df.to_csv("../../data/8-escape_data/4-haplotype_summary_nonsyn/" + virus + "_" + antibody + " " + str(min_value) +"_freq_by_sample.csv")



