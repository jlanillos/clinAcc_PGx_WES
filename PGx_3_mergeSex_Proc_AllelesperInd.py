# This script merges restricted information (per individual gender, NGS panel and country location) with alleles retrieved in the previous script (PGx_2_mergedAlleles_CCP17_SSV6_together.py)
# It creates a dataframe which contains one row per individual (N=5,001) and coluns with all actionable genes filled with their PGx diplotype information whenever they are carriers
# The output file (haplotypes.csv) serves as input to retrieve some statistics, create some SuppTables with alleles per population
# The output file also serves to translate diplotype information final PGx phenotype information per individual

# Important note: In this script, some PGx alleles have undergone manual curation when required (redundant alleles, compound heterozygous cases...) : TPMT, SLCO1B1, G6PD, CYP2B6
import pandas as pd
import numpy as np
import re
GENDER = pd.read_csv('/path/to/dictionary_gender.csv',sep='\t')
proc = pd.read_csv('/path/to/dictionary_proc.csv',sep='\t')
df = proc.copy()
df['gender'] = df['sample'].map(dict(zip(list(GENDER['annonimousID']),list(GENDER['gender'])))) # F=Female;M=Male
panel = pd.read_csv('/path/to/dictionary_panel',sep='\t')
d = dict(zip(list(panel['sample'].astype(str)),list(panel['panel'])))
df['sample'] = df['sample'].astype(str)
df['panel'] = df['sample'].map(d)

# Alleles (SSV6 and CCP17)
alleles = pd.read_csv('/path/to/Alleles.csv',sep='\t')
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()


genot = list(alleles['onlypositivevalues_carrier_gt'])
samples = list(alleles['onlypositivevalues_carrier_ids'])

alleles['aux'] = alleles['allele'].str.replace('_','-') #To avoid misspelling issues later with G6PD, which contains '_' symbol in their allele definitions
alleles['aux'] = alleles['aux'].str.replace(' ','') #To avoid misspelling issues later with G6PD, which contains ', ' symbol in their allele definitions
alleles['aux'] = alleles['aux'].str.replace(' ','') #To avoid misspelling issues later with G6PD, which contains ', ' symbol in their allele definitions
alleles['aux'] = alleles['aux'].str.replace(' ','').str.replace(',',' ')
alleles['ID2'] = alleles['SYMBOL'] + '_' + alleles['aux']

IDsalleles = list(alleles['ID2'])
IDsflat = list()
for i,j in zip(IDsalleles,samples):
    IDsflat.append([i]*len(j.split(',')))
flatIDs =  [item for sublist in IDsflat for item in sublist]

flatsamples = [item for sublist in samples for item in sublist.split(',')]
flatgenot = [item for sublist in genot for item in sublist.split(';')]

cols = ['allele','samples','genotypes']

flatalleles = pd.DataFrame(dict(zip(cols,[flatIDs,flatsamples,flatgenot])))
flatalleles['samples'] = flatalleles['samples'].str.rstrip('_gt')
flatalleles['gender'] = flatalleles['samples'].map(dict(zip(list(df['sample']),list(df['gender']))))
flatalleles = flatalleles.loc[~flatalleles['genotypes'].str.contains('0/0')].copy()
flatalleles['convert_genotypes'] = flatalleles['genotypes'].apply(lambda x: ','.join(list(set(x.split(',')))).replace('1/1,0/1','1').replace('0/1,1/1','1').replace('0/1','1').replace('1/1','2'))

#Useful info to assign phenotypes in G6PD later --> This search gives 0 results since all males will be 1/1 and thus they will get a 2: flatalleles.loc[(flatalleles['allele'].str.contains('G6PD')) & (flatalleles['convert_genotypes'] == '1') & (flatalleles['gender'] == 'H')]

df['alleles'] = ''
df['sample'] = df['sample'].astype(str)
d = dict()

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

d_countalleles = dict()
for s in list(df['sample']):
    l =  list(flatalleles['allele'].loc[flatalleles['samples'] == s].values)
    alleles_aux = list(flatalleles['allele'].loc[flatalleles['samples'] == s].values)
    genot_aux = list(flatalleles['convert_genotypes'].loc[flatalleles['samples'] == s].values)
    d_aux = dict(zip(alleles_aux,genot_aux))

    for i in d_aux.keys():
        if d_aux[i] == '1':
            d_aux[i] = i.replace('_','_*1/')
        elif d_aux[i] == '2':
            d_aux[i] = i + '/' + i.split('_')[1]

    tenta_alleles = ','.join(list(d_aux.values()))
    if ('SLCO1B1_*1/*15' in tenta_alleles) & ('SLCO1B1_*1/*5' in tenta_alleles):
        d_aux = removekey(d_aux, 'SLCO1B1_*5')
    if ('SLCO1B1_*15/*15' in tenta_alleles) & ('SLCO1B1_*5/*5' in tenta_alleles):
        d_aux = removekey(d_aux, 'SLCO1B1_*5')
    if ('SLCO1B1_*1/*15' in tenta_alleles) & ('SLCO1B1_*5/*5' in tenta_alleles):
        d_aux = removekey(d_aux, 'SLCO1B1_*5')
        d_aux = removekey(d_aux, 'SLCO1B1_*15')
        d_aux['SLCO1B1_*5'] = 'SLCO1B1_*5/*15'

    if ('CYP2B6_*1/*34' in tenta_alleles):
        d_aux = removekey(d_aux, 'CYP2B6_*4')
        d_aux = removekey(d_aux, 'CYP2B6_*9')
        d_aux = removekey(d_aux, 'CYP2B6_*6')
        d_aux = removekey(d_aux, 'CYP2B6_*7')
        d_aux = removekey(d_aux, 'CYP2B6_*36')
        if ('CYP2B6_*1/*13' in tenta_alleles):
            d_aux = removekey(d_aux, 'CYP2B6_*8')
        tenta_alleles = ','.join(list(d_aux.values()))
    else:
        if ('CYP2B6_*1/*7' in tenta_alleles):
            d_aux = removekey(d_aux, 'CYP2B6_*4')
            d_aux = removekey(d_aux, 'CYP2B6_*9')
            d_aux = removekey(d_aux, 'CYP2B6_*6')
            if ('CYP2B6_*1/*13' in tenta_alleles):
                d_aux = removekey(d_aux, 'CYP2B6_*8')
            tenta_alleles = ','.join(list(d_aux.values()))
        elif ('CYP2B6_*1/*36' in tenta_alleles):
            d_aux = removekey(d_aux, 'CYP2B6_*4')
            d_aux = removekey(d_aux, 'CYP2B6_*9')
            d_aux = removekey(d_aux, 'CYP2B6_*6')
            d_aux = removekey(d_aux, 'CYP2B6_*22')
            if ('CYP2B6_*1/*13' in tenta_alleles):
                d_aux = removekey(d_aux, 'CYP2B6_*8')
            tenta_alleles = ','.join(list(d_aux.values()))
        elif ('CYP2B6_*1/*8' in tenta_alleles):
            d_aux = removekey(d_aux, 'CYP2B6_*8')
            tenta_alleles = ','.join(list(d_aux.values()))
    if ('CYP2B6_*1/*4' in tenta_alleles) & ('CYP2B6_*1/*6' in tenta_alleles) & ('CYP2B6_*1/*9' in tenta_alleles):
        d_aux = removekey(d_aux, 'CYP2B6_*4')
        d_aux = removekey(d_aux, 'CYP2B6_*9')
        if ('CYP2B6_*1/*13' in tenta_alleles):
            d_aux = removekey(d_aux, 'CYP2B6_*6')
        tenta_alleles = ','.join(list(d_aux.values()))

    if ('CYP2B6_*4/*4' in tenta_alleles) & ('CYP2B6_*6/*6' in tenta_alleles) & ('CYP2B6_*9/*9' in tenta_alleles):
        d_aux = removekey(d_aux, 'CYP2B6_*4')
        d_aux = removekey(d_aux, 'CYP2B6_*9')
        tenta_alleles = ','.join(list(d_aux.values()))

    if ('CYP2B6_*4/*4' in tenta_alleles) & ('CYP2B6_*1/*6' in tenta_alleles) & ('CYP2B6_*1/*9' in tenta_alleles):
        d_aux = removekey(d_aux, 'CYP2B6_*4')
        d_aux = removekey(d_aux, 'CYP2B6_*9')
        d_aux['CYP2B6_*6'] = 'CYP2B6_*4/*6'
        tenta_alleles = ','.join(list(d_aux.values()))

    if ('CYP2B6_*1/*4' in tenta_alleles) & ('CYP2B6_*1/*6' in tenta_alleles) & ('CYP2B6_*9/*9' in tenta_alleles):
        d_aux = removekey(d_aux, 'CYP2B6_*4')
        d_aux = removekey(d_aux, 'CYP2B6_*9')
        d_aux['CYP2B6_*6'] = 'CYP2B6_*9/*6'
        if ('CYP2B6_*1/*13' in tenta_alleles):
            d_aux = removekey(d_aux, 'CYP2B6_*6')
            d_aux['CYP2B6_*9'] = 'CYP2B6_*1/*9'
        tenta_alleles = ','.join(list(d_aux.values()))
    if ('CYP2B6_*1/*16' in tenta_alleles):
        d_aux = removekey(d_aux, 'CYP2B6_*18')
        tenta_alleles = ','.join(list(d_aux.values()))
    if ('CYP2B6_*1/*13' in tenta_alleles)  & ('CYP2B6_*6/*6' in tenta_alleles):
        d_aux['CYP2B6_*6'] = 'CYP2B6_*1/*6'
        tenta_alleles = ','.join(list(d_aux.values()))
    if ('CYP2B6_*1/*26' in tenta_alleles)  & ('CYP2B6_*6/*6' in tenta_alleles):
        d_aux['CYP2B6_*6'] = 'CYP2B6_*1/*6'
    if ('CYP2B6_*1/*20' in tenta_alleles)  & ('CYP2B6_*1/*16' in tenta_alleles)& ('CYP2B6_*1/*6' in tenta_alleles):
        d_aux = removekey(d_aux, 'CYP2B6_*16')
        d_aux = removekey(d_aux, 'CYP2B6_*20')
        d_aux['CYP2B6_*18'] = 'CYP2B6_*1/*18'
        tenta_alleles = ','.join(list(d_aux.values()))
    if ('TPMT_*1/*3A' in tenta_alleles) & ('TPMT_*1/*3B' in tenta_alleles) & ('TPMT_*1/*3C' in tenta_alleles):
        d_aux = removekey(d_aux, 'TPMT_*3B')
        d_aux = removekey(d_aux, 'TPMT_*3C')
    if ('TPMT_*1/*3A' in tenta_alleles) & ('TPMT_*1/*3B' in tenta_alleles) & ('TPMT_*3C/*3C' in tenta_alleles):
        d_aux = removekey(d_aux, 'TPMT_*3B')
        tenta_alleles = ','.join(list(d_aux.values()))
    if ('TPMT_*3A/*3A' in tenta_alleles) & ('TPMT_*3B/*3B' in tenta_alleles) & ('TPMT_*3C/*3C' in tenta_alleles):
        d_aux = removekey(d_aux, 'TPMT_*3B')
        d_aux = removekey(d_aux, 'TPMT_*3C')
        tenta_alleles = ','.join(list(d_aux.values()))
    if ('G6PD_*1/Asahi' in tenta_alleles) & ('G6PD_*1/A-202A-376G' in tenta_alleles):
        d_aux = removekey(d_aux, 'G6PD_Asahi')
        tenta_alleles = ','.join(list(d_aux.values()))
    if ('G6PD_Asahi/Asahi' in tenta_alleles) & ('G6PD_A-202A-376G/A-202A-376G' in tenta_alleles):
        d_aux = removekey(d_aux, 'G6PD_Asahi')
        tenta_alleles = ','.join(list(d_aux.values()))

    d[s] = ','.join(list(d_aux.values()))
    d_countalleles[s] = len(list(d_aux.values()))


df['alleles'] = df['sample'].map(d)
df['N_alleles'] = df['sample'].map(d_countalleles)

GENES = list(set(list(alleles['SYMBOL'])))
for gene in GENES:
    df[gene] = df['alleles'].apply(lambda x: ','.join([j for j in x.split(',') if gene in j]))

df.to_csv('/outdir/onlySNVS_haplotypes_20210107.csv',sep='\t',index = None)


# In this second part, we included indel variation which was manually analysed in parallel to avoid losing positive cases
# The aim is to include the resulting actionable diplotypes from indel variation
# INDELS:

indels = pd.read_csv('/path/to/ActionableIndels_ALLGENES.csv',sep='\t')
sampledict = pd.read_csv('/path/to/dictionary_sample.csv',sep='|', header = None) # Restricted
sampledict = sampledict.rename(columns={0:'sample',1:'NIM'}).copy()
sampledict['sample'] = sampledict['sample'].astype(str)
d = dict(zip(list(sampledict['ORIGINAL']),list(sampledict['sample'])))
def annonymize(d, l):
    m = ','.join([d[i] for i in l])
    return m
indels['samples'] = indels['samples'].apply(lambda x: annonymize(d, x.split(';')))
indels['zigosity'] = indels['zigosity'].apply(lambda x: x.replace('Heteroz','0/1').replace('Homoz','1/1'))

df = df[['sample','from','from_general','gender','panel','alleles','N_alleles']].copy()
df = df.rename(columns={'alleles':'SNV_alleles','N_alleles':'SNV_N_alleles'}).copy()


alleles = pd.read_csv('/path/to/reference_HAPLOTYPES_20201130_hg38_hg19.csv',sep='\t')
alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['SNV_only'] == 'INDEL') & (alleles['actionable'] == 'Yes')]

genot = list(indels['zigosity'])
samples =  list(indels['samples'])
indels['ID2'] = indels['SYMBOL'] + '_' + indels['Allele']
IDsalleles = list(indels['ID2'])
IDsflat = list()
for i,j in zip(IDsalleles,samples):
    IDsflat.append([i]*len(j.split(',')))
flatIDs =  [item for sublist in IDsflat for item in sublist]
flatsamples = [item for sublist in samples for item in sublist.split(',')]
flatgenot = [item for sublist in genot for item in sublist.split(';')]
cols = ['allele','samples','genotypes']
flatalleles = pd.DataFrame(dict(zip(cols,[flatIDs,flatsamples,flatgenot])))
flatalleles['convert_genotypes'] = flatalleles['genotypes'].apply(lambda x: ','.join(list(set(x.split(',')))).replace('1/1,0/1','1').replace('0/1,1/1','1').replace('0/1','1').replace('1/1','2'))

df['alleles'] = ''
df['sample'] = df['sample'].astype(str)
d = dict()
def removekey(d, key):
    r = dict(d)
    del r[key]
    return r
d_countalleles = dict()
for s in list(df['sample']):
    l =  list(flatalleles['allele'].loc[flatalleles['samples'] == s].values)
    alleles_aux = list(flatalleles['allele'].loc[flatalleles['samples'] == s].values)
    genot_aux = list(flatalleles['convert_genotypes'].loc[flatalleles['samples'] == s].values)
    d_aux = dict(zip(alleles_aux,genot_aux))

    for i in d_aux.keys():
        if d_aux[i] == '1':
            d_aux[i] = i.replace('_','_*1/')
        elif d_aux[i] == '2':
            d_aux[i] = i + '/' + i.split('_')[1]
    d[s] = ','.join(list(d_aux.values()))
    d_countalleles[s] = len(list(d_aux.values()))
df['INDELS_alleles'] = df['sample'].map(d)
df['INDELS_N_alleles'] = df['sample'].map(d_countalleles)

df['alleles'] = df['SNV_alleles'].astype(str) + ',' + df['INDELS_alleles'].astype(str)
df['N_alleles'] = df['SNV_N_alleles'] +  df['INDELS_N_alleles']
df.loc[df['N_alleles'] == 0, ['alleles']] = ''
#df['alleles'] = df['alleles'].apply(lambda x: x.replace('nan,','').replace(',nan',''))
df['alleles'] = df['alleles'].apply(lambda x: x.rstrip(',').lstrip(','))




# Include all GENES, those containing Indels and SNVS (that's why I repeat this step of loading "alleles" dataframe) This prevents badly groupping in 20210105_plotStacked...INDELS.py
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))
for gene in GENES:
    df[gene] = df['alleles'].apply(lambda x: ','.join([j for j in x.split(',') if gene in j]))



# Save the final output:
df.to_csv('/path/to/haplotypes_20210107.csv',sep='\t',index = None)
