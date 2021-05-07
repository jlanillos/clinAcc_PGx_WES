# This script will open each .tsv file and extract from the whole set the variants only in the PGx genes of interest. TSV files are annotated with
# relevant biological (variant consequence and impact, etc) and population (MAFs, clinVar, dbSNP) information from various databases. You may convert your annotated
# VCF files into datatable format (e.g., use last metadata line (#CHROM POS...) as a header an remove previous metadata lines)
# Also, you may modify these and subsequet scripts with the right column names and data of your interest

# "genes" will be a dataframe containing the list of gene names (symbols) that will be extracted from the whole data.
# Each gene is written in a new line and the text file will have no header. In this study, we specified the 11 clinically actionable PGx genes:
# CACNA1S,CYP2B6,CYP2C9,CYP4F2,G6PD,DPYD,NUDT15,RYR1,SLCO1B1,TPMT,UG1TA1

# The script outputs per-sample CSV files conaining only variation data of the selected genes that will be considered further

import sys
import re
import argparse
import pandas as pd
import numpy as np
import numbers
import os

parser=argparse.ArgumentParser(description='Discovery_0_PGx_getPGxselection.py')
parser.add_argument('--fileformat', required=True)
parser.add_argument('--searchdir', required=True)
parser.add_argument('--genes', required=True) # Text file containing the gene symbols to be searched from each tsv file with variants of those genes of interest. Each line should contain a gene symbol
parser.add_argument('--outdir', required=True) # Output directory


args=parser.parse_args()
filenameformat=args.fileformat # filenameformat = '.tsv' # This format will find initiavl VCF files to extract sample names
searchdir = args.searchdir #searchpath = '/path/to/tsv_files'
genesfile = args.genes
outdir = args.outdir


# Function to find all files with the same format in te specific folder
def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file


###############################################################################
# Open the file with the genes of interest and save them into a list
genes = pd.read_csv(genesfile, sep='\t', header = None)
genes = list(genes[0])

# Find all '.vcf' files in the specific directory, done by function called 'files' written above
tsv_files = list()
for file in files(searchdir):
    if filenameformat in file:
        tsv_files.append(file)

for file in tsv_files:
    samplename = file.split('_')[0]
    df = pd.read_csv(file, sep='\t')
    ####old approach: dff = df.loc[df['annonimous_GENE'].str.contains('|'.join(genes))].copy() # The columns named "annonimous_GENE" contains the gene symbols of the file
    df['exact_gene'] = df['annonimous_GENE'].apply(lambda x: bool(set(genes).intersection(x.split(','))))
    dff = df.loc[df['exact_gene'] == True].copy()
    dff.to_csv(outdir + '/' + samplename + '.pgx.csv', sep='\t', index = None)
