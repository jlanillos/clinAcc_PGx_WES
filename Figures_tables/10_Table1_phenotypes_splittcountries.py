# This script helps to create TAble 1 (phenotypes per country)
import pandas as pd
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
import numpy as np

# Include all GENES, those containing Indels and SNVS (that's why I repeat this step of loading "alleles" dataframe) This prevents badly groupping in 20210105_plotStacked...INDELS.py
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
#alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))


df = pd.read_csv('/path/to/phenotypes_20210107.csv',sep='\t')
dff = df.copy()

mask = (dff['N_alleles'] != 0) & (dff['Phenotype_G6PD'] != 'G6PD_Normal')
dff_valid = dff[mask]
dff['N_phenotypes'] = 0
dff['Phenotype'] = dff['Phenotype'].apply(lambda x: ','.join([i for i in x.split(',') if i != 'G6PD_Normal']))
dff.loc[mask, 'N_phenotypes'] = dff_valid['Phenotype'].apply(lambda x: len(x.split(',')))

countries =  list(set(dff['from']))
allcountriesmerge = '|'.join(list(set(dff['from'])))
allcountries =  [allcountriesmerge] +  countries
d = dict()
for i in allcountries:
    d[i] = len(dff.loc[dff['from'].str.contains(i)])


cols = [i for i in list(dff.columns) if ('Phenotype_' in i) and not ('_score' in i)]
cols.sort()
count = 0
for gene in cols:
    g = gene.split('_')[1]
    if gene == 'Phenotype_G6PD':
        gf = dff.loc[dff[gene] != 'G6PD_Normal'].groupby(['from',gene])['sample'].count().reset_index().pivot(index =gene, columns='from', values='sample').reset_index().rename(columns={gene: 'Phenotype'})
        dffaux = dff.loc[dff[gene] != 'G6PD_Normal'].groupby([gene])['sample'].count().reset_index()
        daux = dict(zip(list(dffaux[gene]),list(dffaux['sample'])))
        gf[allcountriesmerge] = gf['Phenotype'].map(daux)
    else:
        gf = dff.groupby(['from',gene])['sample'].count().reset_index().pivot(index =gene, columns='from', values='sample').reset_index().rename(columns={gene: 'Phenotype'})
        dffaux = dff.groupby([gene])['sample'].count().reset_index()
        daux = dict(zip(list(dffaux[gene]),list(dffaux['sample'])))
        gf[allcountriesmerge] = gf['Phenotype'].map(daux)

    gf['Gene'] = gene.lstrip('Phenotype_')
    gf['Total'] = 100* (gf[allcountriesmerge]) / (d[allcountriesmerge])
    gf = gf.drop(columns=[allcountriesmerge])
    for country in countries:
        if country not in list(gf.columns):
            gf[country] = [np.nan]*len(gf)
        else:
            gf[country] = 100 *gf[country] / d[country]
    gf = gf[['Gene', 'Phenotype'] + countries + ['Total']].copy()
    if count == 0:
        GF = gf.copy()
        count = 1
    else:
        GF = pd.concat([GF,gf])


GF.to_csv('/path/to/Figures/Table_statsphenotypes_population.csv',sep='\t',index = None)
