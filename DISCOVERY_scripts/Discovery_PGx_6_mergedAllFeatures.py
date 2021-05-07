# This script takes all intermediate files generated before and merges them into a final dataframe: MERGED_prio1_prio2.csv
# This final dataframe will be used to extract final results and plot figures related to the novel variation discovery section

import pandas as pd
import argparse
import os
parser=argparse.ArgumentParser(description='Merge features obtained in parallel computation previously')
parser.add_argument('--searchpath', required=True)
args=parser.parse_args()
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
#######################################
df_concat = pd.read_csv(searchpath + '/' + 'checkpoint_samples.csv',sep='\t')
eachsamplecols = ['genotype', 'zigosity', 'VAF', 'AlleleRatio', 'AlleleCoverage']
for feature in eachsamplecols:
    df = pd.read_csv('checkpoint_' + feature + '.csv',sep='\t')
    df_concat[feature] = ''
    d = dict(zip(list(df['ID']),list(df[feature])))
    df_concat[feature] = df_concat['ID'].map(d)

samplecols = list(df_concat.columns)
varcols = list(set(samplecols) - set(eachsamplecols) - set(['exact_gene']))
df_concat['samples'] = df_concat['samples'].str.lstrip(';')
df_concat['N_samples'] = df_concat['samples'].apply(lambda x: len(x.split(';')))
df_concat['N_homozygous'] = df_concat['zigosity'].apply(lambda x: str(x).split(';').count('Homoz'))
df_concat['N_heterozygous'] = df_concat['zigosity'].apply(lambda x: str(x).split(';').count('Heteroz'))
df_concat['MAF'] = ((2*df_concat['N_homozygous']) + df_concat['N_heterozygous'])/(2*len(FILES)) # Correction needed for X-chromosome gene (e.g., G6PD) variants (include the gender in this calculation)

finalcols = ['ID'] + eachsamplecols + ['N_samples', 'MAF', 'N_homozygous', 'N_heterozygous', 'samples'] + varcols
df_concat = df_concat[finalcols].copy()

df_concat.to_csv(searchpath + '/' + 'MERGED_prio1_prio2.csv',sep='\t',index=None)
df_concat.to_excel(searchpath + '/' + 'MERGED_prio1_prio2.xlsx')
