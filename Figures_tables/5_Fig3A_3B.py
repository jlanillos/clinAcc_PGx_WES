# Script wh helps to plot Figures 3A and 3B
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Include all GENES, those containing Indels and SNVS (that's why I repeat this step of loading "alleles" dataframe) This prevents badly groupping in 20210105_plotStacked...INDELS.py
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
#alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))


dff = pd.read_csv('/path/to/phenotypes_20210107.csv',sep='\t')
mask = (dff['N_alleles'] != 0) & (dff['Phenotype_G6PD'] != 'G6PD_Normal')
dff_valid = dff[mask]
dff['N_phenotypes'] = 0
dff['Phenotype'] = dff['Phenotype'].apply(lambda x: ','.join([i for i in x.split(',') if i != 'G6PD_Normal']))
dff.loc[mask, 'N_phenotypes'] = dff_valid['Phenotype'].apply(lambda x: len(x.split(',')))

gf = dff.groupby('N_phenotypes')['sample'].count()
GF = {'Nr. phenotypes': list(gf.index), 'Count':100*(gf.values / gf.sum()), 'Group':['(N=5001)']*len(gf)}
GF = pd.DataFrame(GF)
tf = GF.iloc[0:4]
d = {'Nr. phenotypes':'[4,7]', 'Count':sum(GF['Count'].iloc[4:]), 'Group':'(N=5001)'}
tf = tf.append(d, ignore_index=True)


bottom = 0
f, ax1 = plt.subplots(figsize=(2,4))
f.set_size_inches(2.7, 4.0)
for i,j, in zip(list(tf['Count'].values), list(tf['Nr. phenotypes'])):
    ax1.bar('N=5001',i,label=j, bottom = bottom, edgecolor = 'black')
    bottom = bottom + i
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[::-1], labels[::-1], loc='center left',bbox_to_anchor=(1.0, 0.5), title='Nr.alleles', fontsize=14,title_fontsize=14) # title = TITLE,
plt.ylabel('%',fontsize=14)
plt.yticks(np.arange(0, 100,10 ))
plt.subplots_adjust(left=0.23, bottom=0.1, right=0.5, top=0.95, wspace=0.14, hspace=0.24)
plt.savefig('/path/to/Figures/Figure_3A_nrphenotypes.png',format = 'png', dpi = 500)
plt.show()

####################################### FIGURE 3B
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Include all GENES, those containing Indels and SNVS (that's why I repeat this step of loading "alleles" dataframe) This prevents badly groupping in 20210105_plotStacked...INDELS.py
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
#alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))

dff = pd.read_csv('/path/to/phenotypes_20210107.csv',sep='\t')
#dff = dff.loc[dff['from'] == 'ESPAÑA']
mask = (dff['N_alleles'] != 0) & (dff['Phenotype_G6PD'] != 'G6PD_Normal')
dff_valid = dff[mask]
dff['N_phenotypes'] = 0
dff['Phenotype'] = dff['Phenotype'].apply(lambda x: ','.join([i for i in x.split(',') if i != 'G6PD_Normal']))
dff.loc[mask, 'N_phenotypes'] = dff_valid['Phenotype'].apply(lambda x: len(x.split(',')))

GENES.sort()
pct_phenot = list()
for gene in GENES:
    pct_phenot.append(100*(dff.groupby('Phenotype_' + gene)['sample'].count().values.sum() / len(dff)))

f, ax1 = plt.subplots(figsize=(6,3.5))
plt.grid(axis='x')
plt.barh(GENES, [100]*len(GENES), align='center', height=.35, color='tab:grey',label='Actionable phenotype')
plt.barh(GENES, pct_phenot, align='center', height=.35, color='tab:red',label='Actionable phenotype',edgecolor = 'k')
plt.xlim([0,100])
plt.xlabel('% population with pharmacogenetic phenotype (n=5001)', fontsize=12)
plt.subplots_adjust(left=0.130, bottom=0.140, right=0.945, top=0.97, wspace=0.14, hspace=0.24)
#plt.savefig('/path/to/Figures/Fig3B.png',format = 'png', dpi = 500)
plt.savefig('Fig3B.png',format = 'png', dpi = 500)

plt.show()









'''### Figure 2A
cols = ['N_alleles','SNV_N_alleles','INDELS_N_alleles']
gf = df.groupby(cols[0])['sample'].count().reset_index()
gf = gf.rename(columns={'sample':cols[0] + '_all'})
dgf = dict(zip(list(df.groupby(cols[1])['sample'].count().index), list(df.groupby(cols[1])['sample'].count().values)))

plt.subplots_adjust(left=0.10, bottom=0.08, right=0.85, top=0.90, wspace=0.14, hspace=0.24)
plt.xticks(rotation=0)
plt.ylim(0,100)
plt.xlabel('')
plt.show()

plt.xticks(rotation=90)
plt.ylim(0,100)
plt.ylabel('%')
plt.show()'''
