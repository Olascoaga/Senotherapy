# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 17:25:57 2021

@author: Samael Olascoaga
@email: olaskuaga@gmail.com
"""

import pandas as pd
import numpy as np
from random import sample, choices
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

reference_mols = 84
tanimoto = 0.7
reference_degree = 0.392
cycles = 1000000

df = pd.read_csv('Similarity_matrix.csv')
df = df.set_index(df['Unnamed: 0'])
mols = df['Unnamed: 0'].tolist()

def Bootstraping(df, mols):
    s = set(mols)
    random_mols = choices(list(s), k = reference_mols)
    random_sample = df.reindex(index=random_mols, columns=random_mols)
    random_sample = random_sample.mask(random_sample < tanimoto, 0)
    G = nx.DiGraph(random_sample.values)
    
    return nx.average_clustering(G)

degree_list = []
for i in range(cycles):
    degree_list.append(Bootstraping(df, mols))

degree_list = np.asarray(degree_list, dtype=float)    
p = ((degree_list >= reference_degree).sum() / cycles) 
print(p)
sns.set_style("white")
sns.despine()
#sns.distplot(degree_list, kde=False, rug=False)
sns.histplot(degree_list, log_scale=False, fill=False, color='k', bins=23)
sns.despine()
plt.ylabel("Frequency")
plt.xlabel("Degree")
plt.title("Bootstraping degree of similarity chemical networks")