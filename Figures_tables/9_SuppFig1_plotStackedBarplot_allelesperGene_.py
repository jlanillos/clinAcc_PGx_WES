# Create Supp. Figure 1A
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('/path/to/haplotypes_20210107.csv',sep='\t')
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))

fig, (ax1, ax2) = plt.subplots(2,1,figsize=(9,7), sharex=True, sharey=True)
# add a big axis, hide frame https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("%", fontsize=12)
plt.ylabel("Gene (Nr.samples)",labelpad=90, fontsize=12)

dff = df.groupby('SNV_N_alleles')[GENES].count().reset_index().copy()
dff = dff.loc[dff['SNV_N_alleles'] != 0].copy()
d = dict()
dtext = dict()
for i in GENES:
    vals = (dff[i].values / np.sum(dff[i].values))*100
    d[i] = vals
    dtext[i] = i + ' (N=' + str(np.sum(dff[i].values)) + ')'
l = list(np.arange(1,np.array(df.groupby('SNV_N_alleles').count().index).max()+1))
d['SNV_N_alleles'] = l
gf = pd.DataFrame(d)
GF = gf.set_index('SNV_N_alleles').T.reset_index().copy()
GF['GENE'] = GF['index'].map(dtext)
GF = GF.drop(columns=['index'])
GF.plot(ax=ax1, x='GENE',kind = 'barh',stacked = True, edgecolor = 'black', mark_right = True).legend(loc='center left',bbox_to_anchor=(1.0, 0.5))

ax1.set_ylabel('')



dff = df.groupby('N_alleles')[GENES].count().reset_index().copy()
dff = dff.loc[dff['N_alleles'] != 0].copy()
d = dict()
dtext = dict()
for i in GENES:
    vals = (dff[i].values / np.sum(dff[i].values))*100
    d[i] = vals
    dtext[i] = i + ' (N=' + str(np.sum(dff[i].values)) + ')'
l = list(np.arange(1,np.array(df.groupby('N_alleles').count().index).max()+1))
d['N_alleles'] = l
gf = pd.DataFrame(d)
GF = gf.set_index('N_alleles').T.reset_index().copy()
GF['GENE'] = GF['index'].map(dtext)
GF = GF.drop(columns=['index'])
GF.plot(ax=ax2, x='GENE',kind = 'barh',stacked = True, edgecolor = 'black', mark_right = True).legend(loc='center left',bbox_to_anchor=(1.0, 0.5))
ax2.set_ylabel('')

plt.xlim(0,100)


plt.subplots_adjust(left=0.205, right=0.910, top=0.955, bottom=0.07, hspace=0.2)
plt.savefig('/path/to/Figures/SuppFig1.png',format = 'png', dpi = 500)
plt.show()


#The fraction of samples with actionable alleles per gene (including indels) found with one or multiple actionable findings
