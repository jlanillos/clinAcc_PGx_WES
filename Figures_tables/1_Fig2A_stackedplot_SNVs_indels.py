# Script to create and plot Fig 2A
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import cm
df = pd.read_csv('/path/to/haplotypes_20210107.csv',sep='\t')
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')

alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))
colormap = cm.get_cmap('Dark2', 8)


### Figure 2A: prepare the data
cols = ['N_alleles','SNV_N_alleles','INDELS_N_alleles']
gf = df.groupby(cols[0])['sample'].count().reset_index()
colors = [colormap(i) for i in list(gf['N_alleles'])]

gf = gf.rename(columns={'sample':cols[0] + '_all'})
dgf = dict(zip(list(df.groupby(cols[1])['sample'].count().index), list(df.groupby(cols[1])['sample'].count().values)))
gf[cols[1]] = gf['N_alleles'].map(dgf)
dgf = dict(zip(list(df.groupby(cols[2])['sample'].count().index), list(df.groupby(cols[2])['sample'].count().values)))
gf[cols[2]] = gf['N_alleles'].map(dgf)
cols = ['N_alleles_all','SNV_N_alleles','INDELS_N_alleles']
for i in cols:
    gf[i] = 100 * (gf[i] / gf[i].sum())

GF = gf.set_index('N_alleles').T.reset_index().rename(columns={'index':'Variants included'})
new_cols = ['SNVs+Indels','SNVs','Indels']
d = dict(zip(cols, new_cols))
GF['Variants included'] = GF['Variants included'].map(d)
# Combine samples with >4 alleles
GF['5-7'] = GF[[5,6,7]].sum(axis=1)
GF = GF.drop(columns=[5,6,7])


# Plot the figure
f, ax1 = plt.subplots(figsize=(5,3))
GF.plot(x='Variants included',ax= ax1, kind = 'bar',stacked = True, color=colors, edgecolor = 'black',mark_right = True, fontsize=16)
handles, labels = ax1.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc='center left',bbox_to_anchor=(1.0, 0.5), title='Nr.alleles', fontsize=10)#,title_fontsize=14) # title = TITLE,
plt.subplots_adjust(left=0.10, bottom=0.08, right=0.85, top=0.95, wspace=0.14, hspace=0.24)
plt.xticks(rotation=0, fontsize=12)
plt.yticks( fontsize=12)
plt.ylim(0,100)
plt.ylabel('%', fontsize=14)
plt.xlabel('')
plt.subplots_adjust(left=0.15, bottom=0.085, right=0.805, top=0.95, wspace=0.14, hspace=0.24)
plt.savefig('/path/to/Figures/Figure_2A_nrAlleles.png',format = 'png', dpi = 500)
plt.show()
