# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 13:14:57 2021

@author: Samael Olascoaga
@email: olaskuaga@gmail.com
"""

import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy.plot import barplot, dotplot
import numpy as np
import seaborn as sns

sns.set_style("whitegrid")

gene_list = pd.read_csv('common.csv', header=None)
glist = gene_list.squeeze().str.strip().tolist()

names = gp.get_library_name() 

enr = gp.enrichr(gene_list= glist,
                 gene_sets=['KEGG_2019_Human'],
                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                 description='KEGG common targets',
                 # no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                )

resultados = enr.results.head(15)
resultados['-log10(FDR)'] = -np.log10(resultados['Adjusted P-value'])
resultados['Genes'] = resultados['Genes'].str.split(';')
resultados['Genes'] = resultados['Genes'].apply(lambda x: len(x))

g = sns.scatterplot(data=resultados, x="-log10(FDR)", y="Term", hue='-log10(FDR)', palette="seismic"
                , size="Genes", sizes=(30, 300), legend=True)

g.legend(loc=6, bbox_to_anchor=(1, 0.5), ncol=1)
plt.ylabel('')
plt.xlabel('-log10(FDR)')
plt.title('KEGG Common targets')

plt.savefig(r'KEGG_common' + '.svg', format='svg', dpi=600, bbox_inches = "tight" )
