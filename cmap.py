# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 17:37:27 2021

@author: Samael Olascoaga
@email: olaskuaga@gmail.com
"""

from cmapPy.pandasGEXpress.parse_gct import parse
import glob
import pandas as pd

gct_files = glob.glob('*.gct') 
mol_names = [] 
for filename in gct_files:
    mol_names.append(filename.split(".gct")[0])

def process(data):
    df = data.data_df.filter(like='24H', axis=1)
    if not len(df.columns) == 1:
        df = df.filter(like='10', axis=1)
        df = df.filter(like='CPC0', axis=1)
    else:
        pass
    return df

col_names = []    
dfs = []
for i in range(0, len(gct_files)):
    data = parse(gct_files[i])
    df = process(data)
    df = pd.DataFrame(df.mean(axis=1))
    if len(df.columns) > 0:
        col_names.append(mol_names[i])
    else:
        pass
    dfs.append(df)

dfs = [df.set_index(data.row_metadata_df['pr_gene_symbol']) for df in dfs] 
df = pd.concat(dfs, axis=1)
df.columns = col_names