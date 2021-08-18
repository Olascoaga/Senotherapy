# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 15:21:49 2021

@author: Samael Olascoaga
@email: olaskuaga@gmail.com
"""

import pandas as pd
import glob
import seaborn as sns

sns.set_context("paper", font_scale=0.5)

files = glob.glob('*.csv') 

dfs = [] # Lista vacia en donde guardaremos los df
for filename in files: # Recorremos los archivos y guardamos cada uno en la lista dfs[]
    data = pd.read_csv(filename,  index_col=[0])
    dfs.append(data)

df = pd.concat(dfs)
DB = df.pop('DB')
copy = df
copy = copy.T
corrMatrix = copy

lut = dict(zip(DB.unique(), ['#0066FF','#FF6600']))
row_colors = DB.map(lut)


g = sns.clustermap(corrMatrix, row_colors=row_colors, 
                         xticklabels=True, yticklabels=True, 
                         cmap='coolwarm', annot=False, fmt=".2f", 
                         dendrogram_ratio=(.1, .2),
                         cbar_pos=(1, .2, .03, .4))

g.ax_heatmap.tick_params(axis='x', rotation=90)
g.ax_heatmap.tick_params(axis='y', rotation=0)
g.fig.savefig(r'similarity' + '.svg', format='svg', dpi=600, bbox_inches="tight")
