import pandas as pd
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt

# Include all GENES, those containing Indels and SNVS (that's why I repeat this step of loading "alleles" dataframe) This prevents badly groupping in 20210105_plotStacked...INDELS.py
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
#alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))
GENES.sort()

df = pd.read_csv('/path/to/phenotypes_20210107.csv',sep='\t')
dff = df.loc[df['from_general'].str.contains('Spain|LATAM')].copy()

N = 5001
N_espana = len(dff.loc[dff['from_general'] == 'Spain'])
N_latam = len(dff.loc[dff['from_general'] == 'LATAM'])
d_N_aux = {'Spain':N_espana, 'LATAM':N_latam}
chi_dict = dict()
chi2score_dict = dict()
alfreq_españa = dict()
alfreq_latam = dict()
alfreq_noalleles_españa = dict()
alfreq_noalleles_latam = dict()

for gene in GENES:

    if gene != 'G6PD':
        contigency= pd.crosstab(dff['from_general'], dff[gene])
        posalleles = list(set([i.split('_')[1].split(',')[0].split('/')[1] for i in list(set(','.join(list(contigency.columns)).split(',')))]))#posalleles = list(set([i.split('_')[1].split(',')[0].split('/')[1] for i in  list(contigency.columns)]))
        for al in posalleles:
            homozall = [i for i in list(contigency.columns) if '/'.join([al,al]) in i]
            heterozall = [i for i in list(contigency.columns) if (al in i) and ~('/'.join([al,al]) in i)] #[i for i in list(contigency.columns) if '/'.join(['*1',al]) in i]
            origin = list((contigency[heterozall].sum(axis=1) + contigency[homozall].sum(axis=1)*2).index)
            vals = list((contigency[heterozall].sum(axis=1) + contigency[homozall].sum(axis=1)*2).values)
            auxdict = dict(zip(origin,vals))
            alfreq_españa[gene + '_' + al] = auxdict['Spain']
            alfreq_latam[gene + '_' + al] = auxdict['LATAM']
            alfreq_noalleles_españa[gene + '_' + al] = 2*N_espana - auxdict['Spain']
            alfreq_noalleles_latam[gene + '_' + al] = 2*N_latam - auxdict['LATAM']
            cont = pd.DataFrame({'from_general':origin, 'N_alleles':vals})
            cont['No_alleles'] = 2*cont['from_general'].map(d_N_aux)
            cont['No_alleles'] = cont['No_alleles'] - cont['N_alleles']
            cont = cont.set_index('from_general')
            chi2, p, dof, ex = chi2_contingency(cont)
            chi2score_dict[gene + '_' + al] = chi2
            chi_dict[gene + '_' + al] = p

    else:
        dff_aux = dff.loc[dff['gender'] == 'M']
        contigency_males= pd.crosstab(dff_aux['from_general'], dff_aux[gene])
        posalleles_males = list(set([i.split('_')[1].split(',')[0].split('/')[1] for i in list(set(','.join(list(contigency_males.columns)).split(',')))]))#posalleles = list(set([i.split('_')[1].split(',')[0].split('/')[1] for i in  list(contigency.columns)]))
        dff_aux = dff.loc[dff['gender'] == 'F']
        contigency_females= pd.crosstab(dff_aux['from_general'], dff_aux[gene])
        posalleles_females = list(set([i.split('_')[1].split(',')[0].split('/')[1] for i in list(set(','.join(list(contigency_females.columns)).split(',')))]))#posalleles = list(set([i.split('_')[1].split(',')[0].split('/')[1] for i in  list(contigency.columns)]))
        g6pdalleles = list(set(posalleles_males + posalleles_females))
        g6pd_aux = dict()
        #g6pd_aux_noalleles = dict()
        for z in g6pdalleles:
            g6pd_aux[z + '_' + 'Spain'] = 0
            g6pd_aux[z + '_' + 'LATAM'] = 0
            #g6pd_aux_noalleles[z + '_' + 'Spain'] = 0
            #g6pd_aux_noalleles[z + '_' + 'LATAM'] = 0
        for al in posalleles_males:
            homozall = [i for i in list(contigency_males.columns) if '/'.join([al,al]) in i]
            heterozall = [i for i in list(contigency_males.columns) if (al in i) and not ('/'.join([al,al]) in i)] #[i for i in list(contigency.columns) if '/'.join(['*1',al]) in i]
            origin = list((contigency_males[heterozall].sum(axis=1) + contigency_males[homozall].sum(axis=1)).index)
            vals = list((contigency_males[heterozall].sum(axis=1) + contigency_males[homozall].sum(axis=1)).values)
            d_aux = dict(zip(origin, vals))
            for j in list(d_aux.keys()):
                g6pd_aux[al + '_' + j] = g6pd_aux[al + '_' + j] + d_aux[j]
                #g6pd_aux_noalleles[al + '_' + j] = g6pd_aux_noalleles[al + '_' + j] + d_aux[j]
        for al in posalleles_females:
            homozall = [i for i in list(contigency_females.columns) if '/'.join([al,al]) in i]
            heterozall = [i for i in list(contigency_females.columns) if (al in i) and not ('/'.join([al,al]) in i)] #[i for i in list(contigency.columns) if '/'.join(['*1',al]) in i]
            origin = list((contigency_females[heterozall].sum(axis=1) + contigency_females[homozall].sum(axis=1)*2).index)
            vals = list((contigency_females[heterozall].sum(axis=1) + contigency_females[homozall].sum(axis=1)*2).values)
            d_aux = dict(zip(origin, vals))
            for j in list(d_aux.keys()):
                g6pd_aux[al + '_' + j] = g6pd_aux[al + '_' + j] + d_aux[j]


        for al in list(set(posalleles_males + posalleles_females)):
            alfreq_españa[gene + '_' + al] = g6pd_aux[al + '_' + 'Spain']
            alfreq_latam[gene + '_' + al] = g6pd_aux[al + '_' + 'LATAM']
            alfreq_noalleles_españa[gene + '_' + al] = 2*len(dff.loc[(dff['gender'] == 'F') & (dff['from_general'] == 'Spain')]) + len(dff.loc[(dff['gender'] == 'M') & (dff['from_general'] == 'Spain')]) - g6pd_aux[al + '_' + 'Spain']
            alfreq_noalleles_latam[gene + '_' + al] = 2*len(dff.loc[(dff['gender'] == 'F') & (dff['from_general'] == 'LATAM')]) + len(dff.loc[(dff_aux['gender'] == 'M') & (dff['from_general'] == 'LATAM')]) - g6pd_aux[al + '_' + 'LATAM']
            origin = ['Spain','LATAM']
            vals = [g6pd_aux[al + '_' + 'Spain'],g6pd_aux[al + '_' + 'LATAM']]
            vals_no_alleles = [alfreq_noalleles_españa[gene + '_' + al],alfreq_noalleles_latam[gene + '_' + al]]
            cont = pd.DataFrame({'from_general':origin, 'N_alleles':vals, 'No_alleles': vals_no_alleles})
            cont = cont.set_index('from_general')
            chi2, p, dof, ex = chi2_contingency(cont)
            chi2score_dict[gene + '_' + al] = chi2
            chi_dict[gene + '_' + al] = p



df_chi = pd.DataFrame({'allele':list(chi_dict.keys()), 'chi2_pval':list(chi_dict.values())})
df_chi['chi2_score'] = df_chi['allele'].map(chi2score_dict)
df_chi['N_alleles_ESPANA'] = df_chi['allele'].map(alfreq_españa)
df_chi['N_alleles_LATAM'] = df_chi['allele'].map(alfreq_latam)
df_chi['N_No_alleles_ESPANA'] = df_chi['allele'].map(alfreq_noalleles_españa)
df_chi['N_No_alleles_LATAM'] = df_chi['allele'].map(alfreq_noalleles_latam)


df_chi['FREQ_ESPANA'] = df_chi['N_alleles_ESPANA'] / (2*N_espana)
df_chi['FREQ_LATAM'] = df_chi['N_alleles_LATAM'] / (2*N_latam)
df_chi_aux = df_chi.copy()
thr = (0.05/len(df_chi))#correction?
df_chi = df_chi.loc[df_chi['chi2_pval'] <= thr]

df_chi.to_csv('Figures/significantAllelesCompSpainANDLatam.csv',sep='\t',index = None) # add CPICs freqs NFE (and AMR?) and save to onlysignificant_CPICfreqs.csv
CPICfreqs = pd.read_csv('/path/to/0_manuscript_scripts/testdata/population_gnomad_CPIC/onlysignificant_CPICfreqs.csv',sep='\t')# G6PD CPIC guidelines do not have population freqs info, and allele A202 is extremely rare in gnomad
df_chi['European(CPIC)'] = df_chi['allele'].map(dict(zip(list(CPICfreqs['allele']),list(CPICfreqs['EUROPEAN_CPICfreq']))))
df_chi['Latino(CPIC)'] = df_chi['allele'].map(dict(zip(list(CPICfreqs['allele']),list(CPICfreqs['LATINO_CPICfreq']))))


genecheck = ''
newxticks = list()
oldticks = list(df_chi['allele'])
for i in oldticks:
    if genecheck != i.split('_')[0]:
        genecheck = i.split('_')[0]
        newxticks.append(i.replace('_',' '))
    else:
        newxticks.append(i.split('_')[1])
d_newxticks = dict(zip(oldticks, newxticks))
df_chi['allele'] = df_chi['allele'].map(d_newxticks)
plt.figure()

plt.scatter(df_chi['allele'], df_chi['FREQ_ESPANA'], c = 'white', s = 20, label = 'SPAIN', marker = '^', edgecolor ="red")
plt.scatter(df_chi['allele'], df_chi['FREQ_LATAM'], c = 'white', s = 20, label = 'LATAM', marker = 's', edgecolor ="blue")
plt.scatter(df_chi['allele'], df_chi['European(CPIC)'], c = 'red', s = 20, label = 'European (CPIC)', marker = '.', edgecolor ="red")
plt.scatter(df_chi['allele'], df_chi['Latino(CPIC)'], c = 'blue', s = 20, label = 'Latino(CPIC)', marker = '.', edgecolor ="blue")

plt.subplots_adjust(left=0.12, bottom=0.30, right=0.95, top=0.92, wspace=0.20, hspace=0.20)
plt.xticks(rotation=75)
plt.legend()
#plt.title('SPAIN vs LATAM allele comparison\n(Chi-square p<=' + "{:.2e}".format(thr) + ')', fontsize=8)
plt.ylabel('Frequencies', rotation=90)
plt.subplots_adjust(left=0.095, bottom=0.305, right=0.978, top=0.975, wspace=0.17, hspace=0.42)
plt.savefig('/path/to/Figures/SuppFigure_2_chi2Alleles_CPIC_Spain_Latam.png',format = 'png', dpi = 1000)
plt.show()
