# Script implemented to plot Figure 3C
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
dff = df.loc[df['from_general'].str.contains('Spain|LATAM')].copy()

mask = (dff['N_alleles'] != 0) & (dff['Phenotype_G6PD'] != 'G6PD_Normal')
dff_valid = dff[mask]
dff['N_phenotypes'] = 0
dff['Phenotype'] = dff['Phenotype'].apply(lambda x: ','.join([i for i in x.split(',') if i != 'G6PD_Normal']))
dff.loc[mask, 'N_phenotypes'] = dff_valid['Phenotype'].apply(lambda x: len(x.split(',')))

N = 5001
N_Spain = len(dff.loc[dff['from_general'] == 'Spain'])
N_latam = len(dff.loc[dff['from_general'] == 'LATAM'])
d_N_from = {'Spain': N_Spain, 'LATAM': N_latam}

chi_dict = dict()
chi2score_dict = dict()
phenofreq_Spain = dict()
phenofreq_latam = dict()
phenofreq_nophenotypes_Spain = dict()
phenofreq_nophenotypes_latam = dict()

dff['Phenotype_G6PD'] = dff['Phenotype_G6PD'].replace('G6PD_Normal', np.nan)

for gene in GENES:

    contigency= pd.crosstab(dff['from_general'], dff['Phenotype_' + gene])
    phenotypes = list(contigency.columns)
    for phe in phenotypes:
        origin = list((contigency[phe]).index)
        vals = list((contigency[phe]).values)
        auxdict = dict(zip(origin,vals))
        phenofreq_Spain[phe] = auxdict['Spain']
        phenofreq_latam[phe] = auxdict['LATAM']
        phenofreq_nophenotypes_Spain[phe] = N_Spain - auxdict['Spain']
        phenofreq_nophenotypes_latam[phe] = N_latam - auxdict['LATAM']
        cont = pd.DataFrame({'from_general':origin, 'N_phenotypes':vals})
        cont['N_No_phenotypes'] = cont['from_general'].map(d_N_from)
        cont['N_No_phenotypes'] = cont['N_No_phenotypes'] - cont['N_phenotypes']
        cont = cont.set_index('from_general')
        chi2, p, dof, ex = chi2_contingency(cont)
        chi2score_dict[phe] = chi2
        chi_dict[phe] = p

df_chi = pd.DataFrame({'phenotype':list(chi_dict.keys()), 'chi2_pval':list(chi_dict.values())})
df_chi['chi2_score'] = df_chi['phenotype'].map(chi2score_dict)
df_chi['N_phenotypes_Spain'] = df_chi['phenotype'].map(phenofreq_Spain)
df_chi['N_phenotypes_LATAM'] = df_chi['phenotype'].map(phenofreq_latam)
df_chi['N_No_phenotypes_Spain'] = df_chi['phenotype'].map(phenofreq_nophenotypes_Spain)
df_chi['N_No_phenotypes_LATAM'] = df_chi['phenotype'].map(phenofreq_nophenotypes_latam)


df_chi['FREQ_Spain'] = df_chi['N_phenotypes_Spain'] / (N_Spain)
df_chi['FREQ_LATAM'] = df_chi['N_phenotypes_LATAM'] / (N_latam)
df_chi_aux = df_chi.copy()
df_chi = df_chi.loc[df_chi['chi2_pval'] <= 0.01] #Any correction?
df_chi = df_chi.sort_values(by='phenotype')

plt.subplots(figsize=(4,4))
plt.scatter(df_chi['phenotype'], df_chi['FREQ_Spain'], c = 'blue', s = 20, label = 'Spain')
plt.scatter(df_chi['phenotype'], df_chi['FREQ_LATAM'], c = 'red', s = 20, label = 'LatAm')

plt.xticks(rotation=90,fontsize=9)
plt.legend(bbox_to_anchor=(1.0, 0.5), loc= 'center left', handletextpad=0.2, handlelength=1)
#plt.title('Spain vs LatAm phenotypes comparison\n(Chi-square p<=.01)', fontsize=10)
plt.ylabel('Frequencies', rotation=90, fontsize=12)
plt.subplots_adjust(left=0.145, bottom=0.275, right=0.785, top=0.975, wspace=0.17, hspace=0.42)

plt.savefig('/path/to/Figures/Figure_3C.png',format = 'png', dpi = 500)
plt.show()
