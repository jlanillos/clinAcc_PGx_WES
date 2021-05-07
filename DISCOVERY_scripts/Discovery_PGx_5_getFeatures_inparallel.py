# This script takes the name of one of the column names with the information related to the genotyping information in each individual
# It creates an intermediate file "checkpoint_{FEATURENAME}.csv" which will be later maerged with the rest of this kind
import pandas as pd
import os
import numpy as np
import re
import argparse
parser=argparse.ArgumentParser(description='Get feature info in parallel')
parser.add_argument('--feature', required=True) # features: 'genotype,zigosity,VAF,AlleleRatio,AlleleCoverage,samples'
parser.add_argument('--searchpath', required=True)

args=parser.parse_args()
feature=str(args.feature)
searchpath=args.searchpath
##############################################################################
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

template_file = searchpath + '/' + 'checkpoint_samples.csv'
df_concat = pd.read_csv(template_file, sep='\t', dtype=str)
df_concat[feature] = ''
df_concat = df_concat[['ID', feature]]
df_concat['ID_aux'] = df_concat['ID']
df_concat.set_index('ID',inplace=True)

for file in FILES:
    sampleid = file.split('.')[0]
    fd = pd.read_csv(file,sep='\t')
    fd['genotype'] = fd['genotype'].str.replace('.','./.')
    if feature == 'VAF':
        fd['VAF'] = fd['VAF'].apply(lambda x: np.round(x, decimals=2))
    fd = fd.astype(str)
    fd['ID'] = fd['locus'].str.replace(':','_') + '_' + fd['REF'] + '_' + fd['genotype'].str.split('/').str[1]
    fd = fd.drop_duplicates(subset = ['ID'], keep = 'first')
    df_concat[feature] =  df_concat.apply(lambda x:  ';'.join([x[feature],fd[feature].loc[fd['ID'] == x['ID_aux']].values[0]]) if x['ID_aux'] in list(fd['ID']) else x[feature],axis=1)

df_concat = df_concat.reset_index()
df_concat[feature] = df_concat[feature].str.lstrip(';')
finalcols = ['ID', feature]
df_concat = df_concat[finalcols].copy()
df_concat.to_csv(searchpath + '/' + 'checkpoint_' + feature + '.csv',sep='\t', index = None)

