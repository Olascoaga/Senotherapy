# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:25:31 2021

@author: olask
"""

import pandas as pd
import glob
import seaborn as sns
import matplotlib.pylab as plt 
from sklearn import preprocessing
from sklearn import decomposition
from sklearn.cluster import KMeans # Our clustering algorithm

sns.set_context("paper", font_scale=2)
sns.set_style("whitegrid")

files = glob.glob('*.csv') 

dfs = [] # Lista vacia en donde guardaremos los df
for filename in files: # Recorremos los archivos y guardamos cada uno en la lista dfs[]
    data = pd.read_csv(filename,  index_col=[0])
    dfs.append(data)

df = pd.concat(dfs)

X = df.iloc[:,0:978].values
y = df.iloc[:,978].astype('category').values

scaler = preprocessing.StandardScaler()
scaler.fit(X)
X_scaled_array = scaler.transform(X)
X_scaled = pd.DataFrame(X_scaled_array)

pca = decomposition.PCA(n_components=2)
pca.fit(X_scaled)
X_fit = pca.transform(X_scaled)

# Put principal components into a data frame so we can plot it.
dfpc = pd.DataFrame(X_fit, columns=['pc1', 'pc2'])
dfpc['class'] = y

model = KMeans(n_clusters=4)
model.fit(dfpc.iloc[:,:2])
labels = model.predict(dfpc.iloc[:,:2])
cluster_labels = model.labels_
dfpc['cluster'] = cluster_labels

f, ax = plt.subplots(figsize=(12, 9))
ax = sns.scatterplot(data=dfpc, x='pc1', y='pc2', hue='cluster', s=250, 
                style='class', legend='full', palette="Set1")

plt.legend(markerscale=3, fontsize=20, loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(r'kmeans_MCF7' + '.svg', format='svg', dpi=600, bbox_inches = "tight" )
