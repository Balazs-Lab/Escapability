import os
import pandas as pd

dfs = []

for snakefile in snakemake.input:
    file = pd.read_csv(snakefile)
    file.drop(columns=file.columns[0], axis=1, inplace=True)
    filename = os.path.basename(snakefile)

    exp = filename.split('-')[0]
    sample = filename.split('-')[1]
    week = f"wk{filename.split('-')[2]}"
    virus = filename.split('-')[3].split('.')[0]
    sample_id = f"{exp}-{sample}-{week}"

    file['exp'] = [exp for _ in range(file.shape[0])]
    file['sample'] = [sample for _ in range(file.shape[0])]
    file['week'] = [week for _ in range(file.shape[0])]
    file['virus'] = [virus for _ in range(file.shape[0])]
    file['sample_id'] = [sample_id for _ in range(file.shape[0])]
    dfs.append(file)
    
df_final = pd.concat(dfs, ignore_index=True)
df_final = df_final.reset_index(drop=True)
df_final.to_csv(snakemake.output[0])




        
