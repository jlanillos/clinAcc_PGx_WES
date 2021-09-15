import sys
import re
import argparse
import pandas as pd
import numpy as np
import numbers
import os

parser=argparse.ArgumentParser(description='PGx_1_merged_to_alleles.py')
parser.add_argument('--fileformat', required=True)
parser.add_argument('--searchpath', required=True)
parser.add_argument('--hapref', required=True)
parser.add_argument('--mergedCSV', required=True)


args=parser.parse_args()
filenameformat=args.fileformat # filenameformat = '.vcf'
searchpath = args.searchpath #searchpath = '/path/to/vcf_files.vcf'
haplotypeRef = args.hapref #haplotypeRef = '/path/to/reference_HAPLOTYPES_20201130_hg38_hg19.csv'
mergedCSV = args.mergedCSV # mergedCSV = 'MERGED_GVCFs.csv'

df = pd.read_csv(mergedCSV, sep='\t')

reference = pd.read_csv(haplotypeRef, sep='\t')
reference['actionable'].loc[(reference['SYMBOL'] == 'CYP4F2') & (reference['allele'] == '*2')] = 'Yes'
reference = reference.loc[(reference['actionable'] == 'Yes')].copy()
# Function to find all files with the same format in te specific folder
def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file
###############################################################################
# Find all '.vcf' files in the specific directory, done by function called 'files' written above: Get samples names to finally create the samples_gt variable
vcf_files = list()
for file in files(searchpath):
    if filenameformat in file:
        vcf_files.append(file)

sample_list = [x.split('.')[0].split('_')[0] for x in vcf_files]
samples_gt = [x + '_gt' for x in sample_list]

header_cols = list(reference.columns)
def multiID(CHR, position, ref, alt):
    position_list = str(position).split(',')
    chr_element = ['chr']*len(list(str(position).split(',')))
    chr_list = [str(CHR)]*len(str(position).split(','))
    chr_final = list(map("".join,list(zip(chr_element,chr_list))))
    ref_list = str(ref).split(',')
    alt_list = str(alt).split(',')
    outlist = list(map("_".join, list(zip(chr_final,position_list,ref_list,alt_list))))
    outlist = ','.join(outlist)
    return(outlist)


reference['multiID'] = reference.apply(lambda x: multiID(x['chr'], x['position_hg19'],x['ref'],x['alt']),axis=1)
CARRIER_IDS_OUT = list()
CARRIER_GT_OUT = list()
WT_IDS_OUT = list()
WT_GT_OUT = list()

for index, row in reference.iterrows():
    positions = row['multiID'].split(',')
    carrier_ids = list()
    carrier_gt = list()
    wildtype_ids = list()
    wildtype_gt = list()
    for j in samples_gt: #[x for x in list(df.columns) if '_gt' in x]
        sample = j.split('_')[0]
        df['searchID'] = df['ID'] + '_' + df[sample + '_REF_VAR'] + '_' + df[sample + '_ALT_VAR']
        genotypes = list(df[j].loc[df['searchID'].isin(positions)].values)
        if len(genotypes) == len(positions):
            carrier_ids.append(j)
            carrier_gt.append(','.join(genotypes))
            wildtype_ids.append('')
            wildtype_gt.append('')
        else:
            position = ['_'.join(x.split('_')[0:2]) for x in row['multiID'].split(',')]
            dff =  df[['ID',j]].loc[df['ID'].isin(position)].set_index('ID')
            genotypes = list(dff[j].loc[position].values)
            #genotypes = list(df[j].loc[df['ID'].isin(position)].values) # This is the old way, not taking control over the order of genotypes
            wildtype_ids.append(j)
            wildtype_gt.append(','.join(genotypes))
            carrier_ids.append('')
            carrier_gt.append('')

    CARRIER_IDS_OUT.append(','.join(carrier_ids))
    CARRIER_GT_OUT.append(';'.join(carrier_gt))
    WT_IDS_OUT.append(','.join(wildtype_ids))
    WT_GT_OUT.append(';'.join(wildtype_gt))

reference['carrier_ids'] = CARRIER_IDS_OUT
reference['carrier_gt'] = CARRIER_GT_OUT
reference['wt_ids'] = WT_IDS_OUT
reference['wt_gt'] = WT_GT_OUT

reference['onlypositivevalues_carrier_ids'] = reference['carrier_ids'].apply(lambda x: ','.join(list(filter(None, x.split(',')))))
reference['onlypositivevalues_carrier_gt'] = reference['carrier_gt'].apply(lambda x: ';'.join(list(filter(None, x.split(';')))))
reference['onlypositivevalues_wt_ids'] = reference['wt_ids'].apply(lambda x: ','.join(list(filter(None, x.split(',')))))
reference['onlypositivevalues_wt_gt'] = reference['wt_gt'].apply(lambda x: ';'.join(list(filter(None, x.split(';')))))

reference['count_carrier_ids'] = reference['onlypositivevalues_carrier_ids'].loc[reference['onlypositivevalues_carrier_ids'] != ''].apply(lambda x: len(x.split(',')))#reference['onlypositivevalues_carrier_ids'].apply(lambda x: len(x.split(',')))
reference['count_wt_ids'] = reference['onlypositivevalues_wt_ids'].loc[reference['onlypositivevalues_wt_ids'] != ''].apply(lambda x: len(x.split(',')))#reference['onlypositivevalues_wt_ids'].apply(lambda x: len(x.split(',')))

cols = ['SYMBOL', 'allele', 'function', 'chr', 'position_hg38', 'position_hg19','positionatNG', 'ref', 'alt', 'amino_acid', 'multiID','N_variants','Buscar Func var','carrier_ids', 'carrier_gt', 'wt_ids','wt_gt','onlypositivevalues_carrier_ids','onlypositivevalues_carrier_gt','onlypositivevalues_wt_ids','onlypositivevalues_wt_gt','count_carrier_ids','count_wt_ids']


# We have identified the samples with different actionable alleles. Save it to a CSV:
reference[cols].to_csv('Alleles' + '.csv',sep='\t',index = None)

# Note this script has been used to create "Alleles.csv" for each panel (SSV6 and SSccp17). They will be used as inputs for the next script (PGx_2_merdAlleles_CCP17_SSV6_together.py)