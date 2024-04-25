import os
import pandas as pd

path = '../data/raw_fastq_files'
metadata = '../data/metadata/experiment_data.csv'
df = pd.read_csv(metadata)
reference_seqs = list(set(df.virus))

def rename(filename):
    experimentID = filename.split('-')[0]

    if any('wk' in identifier for identifier in filename.split('-')):
        week = [identifier.split('wk')[1] for identifier in filename.split('-') if identifier.startswith('wk')][0]
        mouse_num = int([identifier for identifier in filename.split('-') if identifier.startswith('m')][0].split('_')[0][1:])
        mouseID = f"m{mouse_num:03}"
    else:
        week = filename.split('-')[-1][0]
        mouse_num = int([identifier for identifier in filename.split('-') if identifier.startswith('m')][0][1:])
        mouseID = f"m{mouse_num:03}"

    remainder = '_'.join(filename.split('_')[1:])
    virus_ref_seq = df.loc[df.id == f"{experimentID}-{mouseID}", 'virus'].values[0]
    new_name = f"{experimentID}-{mouseID}-{week}-{virus_ref_seq}_{remainder}"

    return new_name

for filename in os.listdir(path):
    if any(seq in filename for seq in reference_seqs):
        continue
    else:
        new_name = rename(filename)
        os.rename(f"{path}/{filename}", f"{path}/{new_name}")
