# -*- coding: utf-8 -*-
"""
Created on Mon May 10 19:03:44 2021

@author: olask
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv('drugbank.csv')

overlap = []
for i in range(0, 1000000):
    set1 = set(df['ID'].sample(n=550, replace=True))
    set2 = set(df['ID'].sample(n=409, replace=True))
    overlap.append(len(set1.intersection(set2)))
    
overlap = np.asarray(overlap, dtype=float)    
p = ((overlap >= 182).sum() / i)

print(p)
sns.set_style("white")
sns.despine()
#sns.distplot(degree_list, kde=False, rug=False)
g = sns.histplot(overlap, log_scale=False, fill=False, color='k', bins=17)
sns.despine()
plt.ylabel("Frequency")
plt.xlabel("Overlap")
#plt.title("")
sns.despine()

fig = g.get_figure()

fig.savefig(r'target_bootstrap' + '.svg', format='svg', dpi=600, bbox_inches="tight")

