# This script opens and grows the "template.csv" file towards a final "template.csv" file which will contain all uniquely identified (ID)
# variants in the whole cohort. Creating this "template.csv" is just a simple and preliminary (and probably unefficient) step to perform the
# real merging operation later on.

import sys
import re
import argparse
import pandas as pd
import argparse
import os

parser=argparse.ArgumentParser(description='Grow the template by appending all possible IDs in the cohort')
parser.add_argument('--file', required=True)
parser.add_argument('--searchpath', required=True)

args=parser.parse_args()
file=args.file # e.g., Sample1.pgx.csv
searchpath=args.searchpath
###############################################################################
template_file = searchpath + '/' + 'template.csv'

df_concat = pd.read_csv(template_file,sep='\t',dtype=str)
samplecols = [i.replace('_merged','') for i in list(df_concat.columns)]
eachsamplecols = ['genotype', 'zigosity', 'VAF', 'AlleleRatio', 'AlleleCoverage','samples']
varcols = list(set(samplecols) - set(eachsamplecols) - set(['exact_gene']))

sampleid = file.split('.')[0]
fd = pd.read_csv(file,sep='\t', dtype=str)
if 'Locus' in list(fd.columns):
    fd = fd.rename(columns={'Locus':'locus'})
    fd.to_csv(file,sep='\t',index = None)
fd['samples'] = [sampleid]*len(fd)
fd['ID'] = fd['locus'].str.replace(':','_') + '_' + fd['REF'] + '_' + fd['genotype'].str.split('/').str[1]
fd['ID_merged'] = fd['locus'].str.replace(':','_') + '_' + fd['REF'] + '_' + fd['genotype'].str.split('/').str[1]
# Loop for annotations should start here
df_concat = df_concat.merge(fd, on = 'ID_merged', how = 'outer')
df_concat = df_concat.drop(samplecols, axis = 1).copy()

df_concat.to_csv(template_file,sep='\t',index = None)
del df_concat
