####################
## SSV6 and CCP17 ##
####################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# x is a vector of coverage values to plot
x = list(np.arange(0,110,10))
# Load reference alleles and extract hg19 "actionable" positions we aim to calculate coverage statistics for each panel (CCP17 and SSV6)
ref = pd.read_csv('reference_HAPLOTYPES_hg38_hg19.csv',sep='\t')
ref['actionable'].loc[(ref['SYMBOL'] == 'CYP4F2') & (ref['allele'] == '*2')] = 'Yes'# This line could be commented to not include CYP4F2*2 in QC
ref['SNV_only'].loc[(ref['SYMBOL'] == 'CYP4F2') & (ref['allele'] == '*2')] = 'SNV'# This line could be commented to not include CYP4F2*2 in QC

actref = ref.loc[(ref['actionable'] == 'Yes')] # & (ref['SNV_only'] == 'SNV')] # & (ref['SYMBOL'] != 'CYP2B6')]
actref['position_hg19'] = actref['position_hg19'].astype(str).apply(lambda x: x.replace('_notfound','').split('_')[0])
posref =[str(x) for x in list(set(list(','.join(actref['position_hg19']).split(','))))]
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
actref['multiID'] = actref.apply(lambda x: multiID(x['chr'], x['position_hg19'],x['ref'],x['alt']),axis=1)


# Get MERGED matrix from SSV6 panel and filter only "actionable" positions in posref
df_ssv6 = pd.read_csv('/path/to/SSV6_panel/QC_MERGED_GVCFs.csv',sep='\t')
df_ssv6_aux = df_ssv6.copy()
df_ssv6 = df_ssv6.loc[df_ssv6['POS'].astype(str).str.contains('|'.join(posref))].copy()
colscov_ssv6 = [i for i in list(df_ssv6.columns) if '_dp' in i]
# Get MERGED matrix from CP17 panel and filter only "actionable" positions in posref
df_ccp17 = pd.read_csv('/path/to/CCP17_panel/QC_MERGED_GVCFs.csv',sep='\t')
df_ccp17_aux = df_ccp17.copy()
df_ccp17 = df_ccp17.loc[df_ccp17['POS'].astype(str).str.contains('|'.join(posref))].copy()
colscov_ccp17 = [i for i in list(df_ccp17.columns) if '_dp' in i]


nrows = int(len(set(list(actref['SYMBOL']))))
ncols = 2
x3 = list()
for i in list(np.arange(0, np.ceil( nrows / 2), 1)):
    x3.append((i,0))
    x3.append((i,1))
d = dict()

genes_sorted = list(set(list(actref['SYMBOL'])))
genes_sorted.sort()

gene_names = ['|'.join(list(set(list(actref['SYMBOL']))))] + genes_sorted #all genes + pergene
d = dict(zip(gene_names, x3))
d['All genes'] = d.pop('|'.join(list(set(list(actref['SYMBOL'])))))


# Store statistical results coverage values, mean, stdev...) into dictionaries
# This part takes a while, consider it for optimization
d2_ssv6 = dict()
d2_ccp17 = dict()
for GENE in gene_names:
    GENEposref =[str(j) for j in list(set(list(','.join(actref['position_hg19'].loc[actref['SYMBOL'].str.contains(GENE)]).split(','))))]
    # Indels
    aux = [str(j) for j in list(set(','.join(actref['multiID'].loc[(actref['SYMBOL'].str.contains(GENE)) & (actref['SNV_only'] == 'INDEL' )]).split(',')))]
    if (aux[0] == '') and (len(aux) == 1):
        GENEposref_indels = []
    else:
        GENEposref_indels = aux
    #SNVs
    aux = [str(j) for j in list(set(','.join(actref['multiID'].loc[(actref['SYMBOL'].str.contains(GENE)) & (actref['SNV_only'] == 'SNV' )]).split(',')))]
    if (aux[0] == '') and (len(aux) == 1):
        GENEposref_snvs = []
    else:
        GENEposref_snvs = aux
    GENE = GENE.replace('|'.join(list(set(list(actref['SYMBOL'])))),'All genes')

    cov_ssv6 = df_ssv6[colscov_ssv6].loc[df_ssv6['POS'].astype(str).str.contains('|'.join(GENEposref))]
    y_ssv6 = list()
    e_ssv6 = list()
    for j in x:
        y_ssv6.append(cov_ssv6.apply(lambda y: sum(map(lambda i: i > j, y)), axis = 1).mean())
        e_ssv6.append(cov_ssv6.apply(lambda y: sum(map(lambda i: i > j, y)), axis = 1).std())
    y2_ssv6 = list(np.array(y_ssv6)/len(cov_ssv6.columns))
    e2_ssv6 = list(np.array(e_ssv6)/len(cov_ssv6.columns))
    d2_ssv6[GENE] = [y2_ssv6, e2_ssv6]
    #CCP17
    cov_ccp17 = df_ccp17[colscov_ccp17].loc[df_ccp17['POS'].astype(str).str.contains('|'.join(GENEposref))]
    y_ccp17 = list()
    e_ccp17 = list()
    for j in x:
        y_ccp17.append(cov_ccp17.apply(lambda y: sum(map(lambda i: i > j, y)), axis = 1).mean())
        e_ccp17.append(cov_ccp17.apply(lambda y: sum(map(lambda i: i > j, y)), axis = 1).std())
    y2_ccp17 = list(np.array(y_ccp17)/len(cov_ccp17.columns))
    e2_ccp17 = list(np.array(e_ccp17)/len(cov_ccp17.columns))
    d2_ccp17[GENE] = [y2_ccp17, e2_ccp17]



# Using the statistical results calculated above and stored in dictionaries (y2_ccp17, e2_ccp17 and y2_ssv6, e2_ssv6) create and plot the figure 1B
fig, axs = plt.subplots(int(np.ceil( nrows / 2)), 2, figsize=(4, 10))
ylim = [0,1.1]

for GENE in gene_names:

    GENEposref =[str(j) for j in list(set(list(','.join(actref['position_hg19'].loc[actref['SYMBOL'].str.contains(GENE)]).split(','))))]
    # Indels
    aux = [str(j) for j in list(set(','.join(actref['multiID'].loc[(actref['SYMBOL'].str.contains(GENE)) & (actref['SNV_only'] == 'INDEL' )]).split(',')))]
    if (aux[0] == '') and (len(aux) == 1):
        GENEposref_indels = []
    else:
        GENEposref_indels = aux
    #SNVs
    aux = [str(j) for j in list(set(','.join(actref['multiID'].loc[(actref['SYMBOL'].str.contains(GENE)) & (actref['SNV_only'] == 'SNV' )]).split(',')))]
    if (aux[0] == '') and (len(aux) == 1):
        GENEposref_snvs = []
    else:
        GENEposref_snvs = aux
    GENE = GENE.replace('|'.join(list(set(list(actref['SYMBOL'])))),'All genes')

    coord_col = int(d[GENE][1])
    coord_row = int(d[GENE][0])


    y2_ssv6 = d2_ssv6[GENE][0]
    e2_ssv6 = d2_ssv6[GENE][1]
    LABEL = 'SSv6 (n=4002)'
    axs[coord_row, coord_col].errorbar(np.array(x)+1, y2_ssv6, yerr=e2_ssv6, fmt='o', label=LABEL, ecolor='lightgray',color = 'blue', ms=3)

    y2_ccp17 = d2_ccp17[GENE][0]
    e2_ccp17 = d2_ccp17[GENE][1]

    LABEL = 'SSccp17 (n=998)'
    axs[coord_row, coord_col].errorbar(np.array(x)-1, y2_ccp17, yerr=e2_ccp17, fmt='o', label = LABEL, ecolor='lightgray', color = 'red', ms=3)
    if (len(GENEposref_snvs) > 1) or (len(GENEposref_snvs) == 0):
        GENEposref_snvs_text = ' SNVs'
    elif len(GENEposref_snvs) == 1:
        GENEposref_snvs_text = ' SNV'
    if (len(GENEposref_indels) > 1) or (len(GENEposref_indels) == 0):
        GENEposref_indels_text = ' indels'

    if len(GENEposref) <= 1:
        LOCUS_word = ' locus'
    else:
        LOCUS_word = ' loci'
    GENE = GENE.replace('All genes', 'All\ genes')
    axs[coord_row, coord_col].set_title( r"$\bf{" + str(GENE) + "}$" + '\n' + '(' + str(len(GENEposref_snvs)) + GENEposref_snvs_text + ', ' + str(len(GENEposref_indels)) + GENEposref_indels_text+ ', ' + str(len(GENEposref)) + LOCUS_word + ')', fontdict={'fontsize':8})
    axs[coord_row, coord_col].tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False, labelsize=6)
    axs[coord_row, coord_col].tick_params(axis = 'y', labelsize=6)
    axs[coord_row, coord_col].set_ylim(ylim)

    if coord_col == 0:
        axs[coord_row,coord_col].set_ylabel('%', fontsize=8, labelpad=1)
    else:
        axs[coord_row,coord_col].yaxis.set_visible(False)
    if coord_row == 5:
        axs[coord_row,coord_col].set_xlabel('Coverage', fontsize=8, labelpad=1)

handles, labels = axs[coord_row, coord_col].get_legend_handles_labels()
fig.legend(handles, labels)
axs[5, 0].tick_params(axis = 'x', which = 'both', bottom = True, labelbottom = True, labelsize=6)
axs[5, 1].tick_params(axis = 'x', which = 'both', bottom = True, labelbottom = True, labelsize=6)
plt.subplots_adjust(left=0.095, bottom=0.060, right=0.98, top=0.910, wspace=0.13, hspace=0.400)
plt.savefig('/path/to/Figures/Fig1B.png', dpi = 400)
plt.show()
