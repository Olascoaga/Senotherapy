# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 18:00:56 2021

@author: olask
"""
import pandas as pd
import glob
import seaborn as sns
import matplotlib.pylab as plt 
from sklearn.decomposition import PCA # Needed for dimension reduction
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans # Our clustering algorithm

sns.set_context("paper", font_scale=1)
files = glob.glob('*.csv') 

dfs = [] # Lista vacia en donde guardaremos los df
for filename in files: # Recorremos los archivos y guardamos cada uno en la lista dfs[]
    data = pd.read_csv(filename,  index_col=[0])
    dfs.append(data)

df = pd.concat(dfs)

X = df.iloc[:,0:978].values
y = df.iloc[:,978].astype('category').values

std_X = StandardScaler().fit_transform(X) # normalizing the data 

pca = PCA(n_components=58)
principalComponents = pca.fit_transform(std_X)

# Plotting the variances for each PC
PC = range(1, pca.n_components_+1)
plt.bar(PC, pca.explained_variance_ratio_, color='blue')
plt.xlabel('Principal Components')
plt.ylabel('Variance %')
plt.xticks(PC)
plt.show()

PCA_components = pd.DataFrame(principalComponents)
inertias = []
# Creating 10 K-Mean models while varying the number of clusters (k)
for k in range(1,10):
    model = KMeans(n_clusters=k)
    # Fit model to samples
    model.fit(PCA_components.iloc[:,:3])
    # Append the inertia to the list of inertias
    inertias.append(model.inertia_)

plt.plot(range(1,10), inertias, '-p', color='blue')
plt.xlabel('number of clusters, k')
plt.ylabel('inertia')
plt.show()