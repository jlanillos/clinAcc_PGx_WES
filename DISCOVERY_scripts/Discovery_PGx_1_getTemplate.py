import pandas as pd
import os
import argparse
parser=argparse.ArgumentParser(description='Get template to start merging individuals files')
parser.add_argument('--searchpath', required=True) # features: 'genotype,zigosity,VAF,AlleleRatio,AlleleCoverage,samples'
args=parser.parse_args()
searchpath=str(args.searchpath)

def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file
###############################################################################
# Find all '.csv' files in the specific directory, done by function called 'files' written above
FILES = list()
#searchpath = '
filenameformat = 'pgx.csv'
for file in files(searchpath):
    if filenameformat in file:
        FILES.append(file)
file = FILES[0]
df_concat = pd.read_csv(file,sep='\t')
sampleid = file.split('.')[0]
df_concat = df_concat.rename(columns={'Locus':'locus'}) # Some samples contain a column called "Locus", with l in capital letters, which produces an error
df_concat['ID'] = df_concat['locus'].str.replace(':','_') + '_' + df_concat['REF'] + '_' + df_concat['genotype'].str.split('/').str[1]
samplecols = list(df_concat.columns)
df_concat[[i + '_merged' for i in samplecols]] = df_concat[samplecols]
df_concat = df_concat.drop(samplecols,axis = 1)


eachsamplecols = ['genotype', 'zigosity', 'VAF', 'AlleleRatio', 'AlleleCoverage','samples']
varcols = list(set(samplecols) - set(eachsamplecols) - set(['exact_gene']))


for s in eachsamplecols:
    #df_concat[s + '_merged'] = df_concat[s + '_merged'].str.split(';').apply(lambda x: ';'.join([i for i in x if i != 'nan']))
    df_concat[s  + '_merged'] = ''


df_concat = df_concat.astype(str)
df_concat.to_csv(searchpath + '/' + 'template.csv',sep='\t',index = None)
