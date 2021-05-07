# This script serves to read through each individual '.pgx.csv' file and extracts common variant information collected from external databases
# It creates an intermediate output file called "checkpoint_samples.csv" which will be merged with other intermediate files generated in the next script
import pandas as pd
import os
import numpy as np
import re
import argparse
import os

parser=argparse.ArgumentParser(description='Annotates SAMPLES column and update Var Columns')
parser.add_argument('--searchpath', required=True)
args=parser.parse_args()
searchpath=str(args.searchpath)

def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file
###############################################################################
FILES = list()
filenameformat = 'pgx.csv'
for file in files(searchpath):
    if filenameformat in file:
        FILES.append(file)

df = pd.read_csv(FILES[0],sep='\t')
# Extract column names from the first individual sample
samplecols = list(df.columns)

# df_concat is the "template.csv" file which will be now filled with annotation data related to the variants (info annotated from databases, no related to the individual's variant)
df_concat = pd.read_csv(searchpath + '/' + 'template.csv',sep='\t', dtype=str)
eachsamplecols = ['genotype', 'zigosity', 'VAF', 'AlleleRatio', 'AlleleCoverage']
varcols = list(set(samplecols) - set(eachsamplecols) - set(['exact_gene']))
del globals()['df']

df_concat.set_index('ID',inplace=True)
df_concat['ID_aux'] = list(df_concat.index)
df_concat['samples'] = ''
df_concat[varcols] = ''

# Extract annotation information from each individual file and update it into the df_concat dataframe
for file in FILES:
    sampleid = file.split('.')[0]
    print(sampleid)
    fd = pd.read_csv(file,sep='\t')
    fd = fd.astype(str)
    fd['ID'] = fd['locus'].str.replace(':','_') + '_' + fd['REF'] + '_' + fd['genotype'].str.split('/').str[1]
    fd = fd.drop_duplicates(subset = ['ID'], keep = 'first')
    df_concat.update(fd[['ID'] +varcols].set_index('ID'))
    df_concat['samples'] =  df_concat.apply(lambda x:  ';'.join([x['samples'],sampleid]) if x['ID_aux'] in list(fd['ID']) else x['samples'],axis=1)
    del globals()['fd']

for s in eachsamplecols:
    df_concat[s] = ''

df_concat = df_concat.reset_index()
df_concat['samples'] = df_concat['samples'].str.lstrip(';')

df_concat['N_samples'] = df_concat['samples'].apply(lambda x: len(x.split(';')))
finalcols = ['ID'] + eachsamplecols + ['N_samples', 'MAF', 'samples'] + varcols
df_concat = df_concat[finalcols].copy()
df_concat.to_csv(searchpath + '/' + 'checkpoint_samples.csv',sep='\t', index = None)
