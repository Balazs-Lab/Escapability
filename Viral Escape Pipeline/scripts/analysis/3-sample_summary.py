import pandas as pd
import csv
import numpy as np


### Analysis 3 - Sample Specific Analysis


print("Coverage Per Sample Analysis")
### Quality Check - Average Coverage per Sample
coverage_df = pd.read_csv("../../data/7-merged_output/coverage_merged.csv")
coverage_mean = coverage_df.groupby(['sample_id'])["COVERAGE"].mean().reset_index()
coverage_mean.columns = ["sample_id","mean_coverage"]
coverage_std = coverage_df.groupby(['sample_id'])["COVERAGE"].std().reset_index()
coverage_std.columns = ["sample_id","std_coverage"]
out_df = pd.merge(coverage_mean,coverage_std, how='outer',on='sample_id')
out_df.to_csv("../../data/8-escape_data/6-coverage_summary/coverage_summary.csv")



### Analysis 3 - Site specific mutations by sample

def mutation_by_sample(file_id, virus):
    """
    Input: variant frequency file

    Output: HXB2 Aligned Site Specifc Mutation Frequencies (by AA change)
    also reports the coerage at each position
    """

    #read in coverage and variant data
    coverage_data = pd.read_csv("../../data/6-coverage/"+file_id+"-"+virus+".csv")
    variant_data = pd.read_csv("../../data/5-calls/"+file_id+"-"+virus+".csv")

    #merjge coverage data onto variant dataframe
    escape_df = variant_data[["POS_AA","REF_AA","ALT_AA","ALT_FREQ"]].merge(coverage_data[["POS_AA","COVERAGE"]],how='outer',on='POS_AA').fillna(0)
    escape_df.columns = ["aa_pos","REF_AA","ALT_AA","ALT_FREQ","COVERAGE"]

#    pivot_columns = ["aa_pos","COVERAGE","*","A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    pivot_columns = ["aa_pos","COVERAGE","*","A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","DEL","*i","Ai","Ri","Ni","Di","Ci","Ei","Qi","Gi","Hi","Ii","Li","Ki","Mi","Fi","Pi","Si","Ti","Wi","Yi","Vi"]

    pivot_df = escape_df.pivot(index=["aa_pos","COVERAGE"],columns="ALT_AA",values="ALT_FREQ").fillna(0).reset_index().sort_values("aa_pos").reindex(columns=pivot_columns).fillna(0)


    # align to HXB2 and ensure full range of samples
    alignment = pd.read_excel("../../data/static/HXB2 Alignment.xlsx",sheet_name=virus)
    position = alignment[[virus + " Numbering", "HXB2 Numbering",virus+"_Env","HXB2_Env"]]
    position.columns = ["aa_pos","HXB2_pos",virus+"_aa","HXB2_aa"]
    summary_df = position.merge(pivot_df,how='left',on="aa_pos").fillna(0)
    #summary_df = pivot_df.merge(position,how='right',on="aa_pos").fillna(0)

    summary_df["sum"] = summary_df[summary_df.columns[5:]].sum(axis=1)
    summary_df["WT"] = 1-summary_df["sum"]


    summary_df.to_csv("../../data/escape_data/5-sample_mutations/" + file_id + "_" + virus + "_mut_freq.csv")

# run on all the files var_freq
var_freq_df = pd.read_csv("../../data/7-merged_output/var_freq_merged.csv")
exp_df = pd.read_csv("../../data/static/experiment_data.csv")

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
        print(file_id,virus)
        mutation_by_sample(file_id,virus)


CD91-m737
