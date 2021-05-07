# This script just renames some columns in the "template.csv" file

import pandas as pd
import argparse

parser=argparse.ArgumentParser(description='modify the columns (remove _merged= from the final template file')
parser.add_argument('--searchpath', required=True)

args=parser.parse_args()
searchpath=args.searchpath
###############################################################################
template_file = searchpath + '/' + 'template.csv'
######################################################
df_concat = pd.read_csv(template_file,sep='\t')
df_concat  = df_concat.rename(columns=dict(zip([i for i in list(df_concat.columns)],[i.replace('_merged','') for i in list(df_concat.columns)])))
df_concat.to_csv(template_file,sep='\t',index = None)
