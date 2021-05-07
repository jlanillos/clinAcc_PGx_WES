# After running PGx_1_merged_to_alleles.py script for both merged VCFs with the samples of the two different panels (SSV6 and CCP17), load the output CSV files and merged them
import pandas as pd
df_ssv6 = pd.read_csv('/path/to/SSV6_panel/Alleles.csv',sep='\t')
df_ccp17 = pd.read_csv('/path/to/CCP17_panel/Alleles.csv',sep='\t')
# df_concat will be the dataframe to merge both panels. Copy SSV6 df as df_concat
df_concat = df_ssv6.copy()
df_ssv6['ID'] = df_ssv6['SYMBOL'] + '_' + df_ssv6['allele'] + '_' +  df_ssv6['function']
df_ccp17['ID'] = df_ccp17['SYMBOL'] + '_' + df_ccp17['allele'] + '_' +  df_ccp17['function']
df_concat['ID'] = df_concat['SYMBOL'] + '_' + df_concat['allele'] + '_' +  df_concat['function']

# Create two new columns (['onlypositivevalues_carrier_ids','onlypositivevalues_wt_ids']) to write down carrier sample names and wt sample names for each allele and wipe out thest
features = ['onlypositivevalues_carrier_ids','onlypositivevalues_wt_ids']
for feature in features:
    df_ccp17[feature] = df_ccp17[feature].astype(str)
    df_concat[feature] = df_concat[feature].astype(str)
    df_concat[feature] = df_concat.apply(lambda x:  ','.join([x[feature],df_ccp17[feature].loc[df_ccp17['ID'] == x['ID']].values[0]]) if x['ID'] in list(df_ccp17['ID']) else x[feature],axis=1)
    df_concat[feature] = df_concat[feature].str.replace('nan,nan','')
    df_concat[feature] = df_concat[feature].str.replace('nan,','')
    df_concat[feature] = df_concat[feature].str.replace(',nan','')

# Create two new columns (['onlypositivevalues_carrier_ids','onlypositivevalues_wt_ids']) to write down carrier samples genotypes (gt) and wt samples genotypes (gt) for each allele and wipe out thest

features = ['onlypositivevalues_carrier_gt','onlypositivevalues_wt_gt']

for feature in features:
    df_ccp17[feature] = df_ccp17[feature].astype(str)
    df_concat[feature] = df_concat[feature].astype(str)
    df_concat[feature] = df_concat.apply(lambda x:  ';'.join([x[feature],df_ccp17[feature].loc[df_ccp17['ID'] == x['ID']].values[0]]) if x['ID'] in list(df_ccp17['ID']) else x[feature],axis=1)
    df_concat[feature] = df_concat[feature].str.replace('nan;nan','')
    df_concat[feature] = df_concat[feature].str.replace('nan;','')
    df_concat[feature] = df_concat[feature].str.replace(';nan','')

# Count of the newly added columns
df_concat['count_carrier_ids'] = df_concat['onlypositivevalues_carrier_ids'].loc[df_concat['onlypositivevalues_carrier_ids'] != ''].apply(lambda x: len(x.split(',')))
df_concat['count_wt_ids'] = df_concat['onlypositivevalues_wt_ids'].loc[df_concat['onlypositivevalues_wt_ids'] != ''].apply(lambda x: len(x.split(',')))

# "panel_id" is a new column with the panel of each sample carrying each allele (eg: Allele X carriers: sample1, sample2, sample3 --> panel_id: SSV6,SSV6,CCP17)
feature = 'panel_id'
df_concat[feature] = df_concat['count_carrier_ids'].apply(lambda x: ','.join(['SSV6']*int(str(x).split('.')[0].replace('nan','0'))))
df_ccp17[feature] = df_ccp17['count_carrier_ids'].apply(lambda x: ','.join(['CCP17']*int(str(x).split('.')[0].replace('nan','0'))))
df_concat[feature] =  df_concat.apply(lambda x:  ','.join([x[feature],df_ccp17[feature].loc[df_ccp17['ID'] == x['ID']].values[0]]) if x['ID'] in list(df_ccp17['ID']) else x[feature],axis=1)
df_concat[feature] = df_concat[feature].str.lstrip(',')

features = ['carrier_ids','carrier_gt','wt_ids','wt_gt']
df_concat = df_concat.drop(features,axis=1).copy()

ref = pd.read_csv('/path/to/reference_HAPLOTYPES_hg38_hg19.csv',sep='\t')
ref['ID2'] = ref['SYMBOL'] + '_' + ref['allele']
df_concat['ID2'] = df_concat['SYMBOL'] + '_' + df_concat['allele']
df_concat['actionable'] = df_concat['ID2'].map(dict(zip(list(ref['ID2']),list(ref['actionable']))))
df_concat['old_function'] = df_concat['function']
df_concat['function'] = df_concat['ID2'].map(dict(zip(list(ref['ID2']),list(ref['function']))))

# Save the output into a newly created dataframe which merge the allele data from both panels
df_concat.to_csv('Alleles_20201228.csv',sep='\t',index = None)
