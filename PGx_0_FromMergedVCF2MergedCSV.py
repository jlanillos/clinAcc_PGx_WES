import sys
import re
import argparse
import pandas as pd
import numpy as np
parser=argparse.ArgumentParser(description='Convert a merged VCF file with relevant PGx loci into CSV format')
parser.add_argument('--vcf', required=True) # 'MERGED.vcf'

args=parser.parse_args()
filename=args.vcf
infile = open(filename, 'r')

for l in infile:

#metadata
    if l.startswith('##'): continue     #other '##' lines unused
    if l.startswith('#CHROM'): #'#CHROM' is the header
        HEADER=l.strip('#').split()
        SAMPLES=HEADER[9:] # The sample name is on the 9th (and next) field of the VCF file header.
            # Create the Output HEADER (Filter: it's the one established after merging VCF files with bcftools merge (default); DP_allsamples (sum of the DP of all samples containing that variant (bftools merge -i DP:sum)
        OUTPUT_HEADER = ['ID','CHROM','POS','REF','ALT'] # check the INFO header of your VCF for more ouput options
        SAMPLES_HEADER = []
        for sample in SAMPLES: SAMPLES_HEADER = SAMPLES_HEADER + ['_'.join([sample,'gt'])] # + ['_'.join([sample,'ad'])]+ ['_'.join([sample,'dp'])] + ['_'.join([sample,'af'])]
        AF_samplenames = []
        for sample in SAMPLES: AF_samplenames = AF_samplenames + ['_'.join([sample,'af'])] # SAMPLE_NAMES_af to be used later during filtering
        OUTPUT_HEADER = OUTPUT_HEADER + SAMPLES_HEADER
        OUTPUT= []
        continue
#variants
    # s will contain all fields for this line; Store all fields of a variant in a dictionary
    s=l.strip().split('\t')
    s=dict(zip(HEADER, s))

# Extract variant fields
    ID = '_'.join([s['CHROM'],s['POS']]) # Variant ID containing the chrom, pos

# Extract information per sample
    SAMPLE_INFO_fields = []
    OUTPUT_general = [ID, s['CHROM'],s['POS'],s['REF'],s['ALT']]
    for SAMPLE in SAMPLES:
        SAMPLE_INFO = dict(zip(s['FORMAT'].split(':'),s[SAMPLE].split(':'))) # Ex: {'GT': '0/1', 'AD': '33,2', 'AF': '0.081', ... , 'SA_MAP_AF': '0.00,0.061,0.057', 'SA_POST_PROB': '0.011,0.019,0.970'}
        SAMPLE_INFO_fields = SAMPLE_INFO_fields + [SAMPLE_INFO['GT']]#, SAMPLE_INFO['AD'],SAMPLE_INFO['DP'], SAMPLE_INFO['AF']] # SAMPLE_INFO['AD'][0], SAMPLE_INFO['AD'][1] are ref and alt allelic depths, respectively

    OUTPUT.append(OUTPUT_general + SAMPLE_INFO_fields)

# OUTPUT list stored into pandas df
df = pd.DataFrame(OUTPUT,columns = OUTPUT_HEADER)
df['ID'] = df['ID'] + '_' + df['REF'] + '_' + df['ALT']
df = df.loc[df['ALT'] != '<NON_REF>'] # Haplotypecaller in BP_RESOLUTION mode releases the '<NON_REF>' value if ALT variant is Wild-type (discard them)
df = df[['ID'] + SAMPLES_HEADER] # Select columns of interest and sort them

# Fill the gaps not called as 0/0
df = df.apply(lambda x: x.replace('./.','0/0'),axis = 0).groupby('ID').agg(','.join).reset_index()
df = df.applymap(lambda x: ','.join([i for i in list(set(x.split(','))) if i != '0/0'])).apply(lambda z: z.replace('','0/0'))

df.to_csv('MERGED_GVCFs.csv',sep='\t', index = None)

