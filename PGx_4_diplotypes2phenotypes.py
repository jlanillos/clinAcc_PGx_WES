# This script would serve to convert per individual haplotype/diplotype PGx information (from PGx_3_mergeSex_Proc_Allelesperind.py) into the actionable phenotypes
import pandas as pd

df = pd.read_csv('/path/to/haplotypes_20210107.csv',sep='\t')
l = [str(i).split(',') for i in list(df['alleles'])]
flat = [item for sublist in l for item in sublist]
auxflat = list(set(list(flat)))
auxflat.remove('nan')
hapo = pd.DataFrame({'findings':auxflat})

df = pd.read_csv('/path/to/Haplo2Pheno_convertion_dictionary.csv',sep='\t')
hapo['Phenotype'] = hapo['findings'].map(dict(zip(list(df['findings']),list(df['Phenotype']))))
hapo['Activity_score'] = hapo['findings'].map(dict(zip(list(df['findings']),list(df['Activity_score']))))
hapo.to_csv('/path/to/Haplo2Pheno_convertion_dictionary.csv',sep='\t',index = None)
# Save hapo df and perform a manual annotation of a second column in hapo to create a dictionary with the phenotypes

# Load the manually curated Hapo2Pheno convertion dictionary again (hapo) as well as the haplotypes (df)
df = pd.read_csv('/path/to/haplotypes_20210107.csv',sep='\t')
df['alleles'] = df['alleles'].astype(str)
hapo = pd.read_csv('/path/to/Haplo2Pheno_convertion_dictionary.csv',sep='\t')
hapo['Phenotype_score'] = hapo['Phenotype'] + '_' + hapo['Activity_score'].astype(str)
hapo['Phenotype_score'] = hapo['Phenotype_score'].apply(lambda x: x.rstrip('_.').rstrip('_G6PD'))
hapodict = dict(zip(list(hapo['findings']),list(hapo['Phenotype'])))
hapodict['nan'] = ''
def convert2pheno(hapodict, alleles):
    return ','.join([j.split('_')[0] + '_' + hapodict[j] for j in alleles.split(',')])

df['Phenotype'] = df['alleles'].apply(lambda x: convert2pheno(hapodict, x)).replace('nan','')
haposcoredict = dict(zip(list(hapo['findings']),list(hapo['Phenotype_score'])))
haposcoredict['nan'] = ''
df['Phenotype_score']  = df['alleles'].apply(lambda x: convert2pheno(haposcoredict, x)).replace('nan','')

# Include all GENES, those containing Indels and SNVS (that's why I repeat this step of loading "alleles" dataframe) This prevents badly groupping in 20210105_plotStacked...INDELS.py
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))
cols = list()
scorecols = list()
for gene in GENES:
    df['Phenotype_' + gene] = df['Phenotype_score'].apply(lambda x: ','.join([j for j in x.split(',') if gene in j]))
    cols.append('Phenotype_' + gene)

#Changing specific cases:
# Female with two deficient G6PD alleles:
df['Phenotype_' + 'G6PD'].loc[(df['G6PD'].astype(str).str.contains('Seattle')) & (df['G6PD'].str.contains(',')) & (df['Phenotype_' + 'G6PD'] == 'G6PD_Normal,G6PD_Normal')] = 'G6PD_Deficient'

# CYP2B6
# CYP2B6_tosolve = list(set(list(df['CYP2B6'].loc[(df['Phenotype_CYP2B6'].str.contains(','))].values)))
CYP2B6_tosolve = dict(zip(list(df['CYP2B6'].loc[(df['Phenotype_CYP2B6'].str.contains(','))].values),df['Phenotype_CYP2B6'].loc[(df['Phenotype_CYP2B6'].str.contains(','))].values))
print(CYP2B6_tosolve)
# Manually corrected: CYP2B6_tosolve = {'CYP2B6_*1/*7,CYP2B6_*1/*8': 'CYP2B6_IM,CYP2B6_IM',
# 'CYP2B6_*1/*6,CYP2B6_*1/*18': 'CYP2B6_IM,CYP2B6_IM',
# 'CYP2B6_*1/*18,CYP2B6_*1/*20': 'CYP2B6_IM,CYP2B6_IM',
# 'CYP2B6_*1/*9,CYP2B6_*1/*22': 'CYP2B6_IM,CYP2B6_RM',
# 'CYP2B6_*1/*6,CYP2B6_*1/*12': 'CYP2B6_IM,CYP2B6_IM',
# 'CYP2B6_*1/*4,CYP2B6_*1/*22': 'CYP2B6_RM,CYP2B6_RM',
# 'CYP2B6_*1/*8,CYP2B6_*1/*36': 'CYP2B6_IM,CYP2B6_IM'}
# The final convertion will be:{'CYP2B6_IM,CYP2B6_IM': 'CYP2B6_PM', 'CYP2B6_IM,CYP2B6_RM': 'CYP2B6_IM', 'CYP2B6_RM,CYP2B6_RM': 'CYP2B6_URM'}
df['Phenotype_' + 'CYP2B6'] = df['Phenotype_' + 'CYP2B6'].astype(str).apply(lambda x: x.replace('CYP2B6_IM,CYP2B6_IM', 'CYP2B6_PM').replace('CYP2B6_IM,CYP2B6_RM', 'CYP2B6_IM').replace('CYP2B6_RM,CYP2B6_RM', 'CYP2B6_URM'))

#CYP2C9
CYP2C9_tosolve = dict(zip(list(df['CYP2C9'].loc[(df['Phenotype_CYP2C9'].str.contains(','))].values),df['Phenotype_CYP2C9'].loc[(df['Phenotype_CYP2C9'].str.contains(','))].values))
#CYP2C9_tosolve= {'CYP2C9_*1/*2,CYP2C9_*1/*12': 'CYP2C9_IM_1.5,CYP2C9_IM_1.5',
# 'CYP2C9_*1/*2,CYP2C9_*1/*3': 'CYP2C9_IM_1.5,CYP2C9_PM_1',
# 'CYP2C9_*1/*2,CYP2C9_*1/*11': 'CYP2C9_IM_1.5,CYP2C9_IM_1.5',
# 'CYP2C9_*1/*2,CYP2C9_*1/*14': 'CYP2C9_IM_1.5,CYP2C9_IM_1.5',
# 'CYP2C9_*1/*3,CYP2C9_*1/*11': 'CYP2C9_PM_1,CYP2C9_IM_1.5',
# 'CYP2C9_*1/*2,CYP2C9_*1/*8': 'CYP2C9_IM_1.5,CYP2C9_IM_1.5',
# 'CYP2C9_*1/*3,CYP2C9_*1/*8': 'CYP2C9_PM_1,CYP2C9_IM_1.5',
# 'CYP2C9_*1/*2,CYP2C9_*1/*6': 'CYP2C9_IM_1.5,CYP2C9_PM_1',
# 'CYP2C9_*1/*3,CYP2C9_*1/*12': 'CYP2C9_PM_1,CYP2C9_IM_1.5'}
# The final convertion will be: {'CYP2C9_IM_1.5,CYP2C9_IM_1.5':'CYP2C9_IM_1','CYP2C9_IM_1.5,CYP2C9_PM_1':'CYP2C9_PM_0.5','CYP2C9_PM_1,CYP2C9_IM_1.5':'CYP2C9_PM_0.5'}
df['Phenotype_' + 'CYP2C9'] = df['Phenotype_' + 'CYP2C9'].astype(str).apply(lambda x: x.replace('CYP2C9_IM_1.5,CYP2C9_IM_1.5','CYP2C9_IM_1').replace('CYP2C9_IM_1.5,CYP2C9_PM_1','CYP2C9_PM_0.5').replace('CYP2C9_PM_1,CYP2C9_IM_1.5','CYP2C9_PM_0.5'))

#DPYD
DPYD_tosolve = dict(zip(list(df['DPYD'].loc[(df['Phenotype_DPYD'].str.contains(','))].values),df['Phenotype_DPYD'].loc[(df['Phenotype_DPYD'].str.contains(','))].values))
#DPYD_tosolve = {'DPYD_*1/*2A,DPYD_*1/HapB3': 'DPYD_IM_1,DPYD_IM_1.5'}
# The final convertion will be: {'DPYD_IM_1,DPYD_IM_1.5': 'DPYD_PM_0.5'}
df['Phenotype_' + 'DPYD'] = df['Phenotype_' + 'DPYD'].astype(str).apply(lambda x: x.replace('DPYD_IM_1,DPYD_IM_1.5', 'DPYD_PM_0.5'))

#TPMT
TPMT_tosolve = dict(zip(list(df['TPMT'].loc[(df['Phenotype_TPMT'].str.contains(','))].values),df['Phenotype_TPMT'].loc[(df['Phenotype_TPMT'].str.contains(','))].values))
#TPMT_tosolve = {'TPMT_*1/*3A,TPMT_*1/*2': 'TPMT_IM,TPMT_IM', 'TPMT_*1/*3A,TPMT_*3C/*3C': 'TPMT_IM,TPMT_PM'}
# The final convertion will be: {'TPMT_IM,TPMT_IM': 'TPMT_PM', 'TPMT_IM,TPMT_PM': 'TPMT_PM'}
df['Phenotype_' + 'TPMT'] = df['Phenotype_' + 'TPMT'].astype(str).apply(lambda x: x.replace('TPMT_IM,TPMT_IM', 'TPMT_PM').replace('TPMT_IM,TPMT_PM', 'TPMT_PM'))

UGT1A1_tosolve = dict(zip(list(df['UGT1A1'].loc[(df['Phenotype_UGT1A1'].str.contains(','))].values),df['Phenotype_UGT1A1'].loc[(df['Phenotype_UGT1A1'].str.contains(','))].values))
#UGT1A1_tosolve = {'UGT1A1_*1/*6,UGT1A1_*1/*28': 'UGT1A1_IM,UGT1A1_IM', 'UGT1A1_*1/*27,UGT1A1_*1/*28': 'UGT1A1_IM,UGT1A1_IM'}
# The final convertion will be: {'UGT1A1_IM,UGT1A1_IM':'UGT1A1_IM'}
df['Phenotype_' + 'UGT1A1'] = df['Phenotype_' + 'UGT1A1'].astype(str).apply(lambda x: x.replace('UGT1A1_IM,UGT1A1_IM', 'UGT1A1_IM'))


for pheno in cols:
    gene = pheno.split('_')[1]
    print(gene)
    print(list(set(list(df.loc[df[pheno].str.contains(',')][gene]))))
    print('\n')


# The output file is a CSV with per individual PGx phenotype information

df.to_csv('/path/to/phenotypes_20210107.csv',sep='\t',index = None)
