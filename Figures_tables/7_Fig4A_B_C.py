# Script to plot Figures 4 (A, B and C)
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# Prepare the dataframe containing all variation data. MERGED_prio1_prio2.csv is a dataframe with all germline variation found in actionable genes (known and novel)

df = pd.read_csv('/path/to/MERGED_prio1_prio2.csv',sep='\t') # '/mnt/64716603-5b56-4f9a-b195-c11560647a3a/Projects/PHARMACOGENETICS/PGx_project/CNIO_jLanillos/Tier1_PharmGKB/samples'
#Filtering the "consequence" column for those variants of interest. I have previously checked all possible consequences annotated and chosen the following ones:
df = df.loc[df['annonimous_ANNOTATION'].str.contains('missense|frameshift|stop|start_lost|splice')] # len(df) = 3448
mask = (df['annonimous_GENE'].str.contains(','))
df_valid = df[mask]
df['SYMBOL'] = df['annonimous_GENE']
# Some variants affect overlapping genes, and I just want to get the info related to the genes of interest. Create an auxiliary "'aux'" column which contains the index location of the gene of interest. E.g. UGT1A8,UGT1A4,UGT1A1
# aux = 2, to retrieve latter the gene UGT1A1 from that list (which is in "annonimous_GENE" col)
genesdf = pd.read_csv('/path/to/bioMart_transcripts_length.csv',sep='\t')
genes = list(genesdf['Gene name'])
df.loc[mask, 'aux'] = df_valid.apply(lambda x: str([i for i, j in enumerate(x['annonimous_GENE'].split(',')) if j == [z for z in genes if z in x['annonimous_GENE']][0]]).replace('[','').replace(']',''), axis=1)
mask = (df['annonimous_GENE'].str.contains(','))
df_valid = df[mask]
df.loc[mask, 'SYMBOL'] = df_valid.apply(lambda x: x['annonimous_GENE'].split(',')[int(x['aux'])],axis = 1)
# Since we know the 'aux' column, we can apply the same principle to columns with multiple annotations and get the right term into "ANNOTATION" column
df['ANNOTATION'] = df['annonimous_ANNOTATION']
df.loc[mask, 'ANNOTATION'] = df_valid.apply(lambda x: x['annonimous_ANNOTATION'].split(',')[int(x['aux'])],axis = 1)
#Filter by consequence again on the newly created "ANNOTATION" column. That column may contain a consequence we did not want
df = df.loc[df['ANNOTATION'].str.contains('missense|frameshift|stop|start_lost|splice')] # len(df) = 3387
df = df.loc[~(df['ANNOTATION'].str.contains('synonymous SNV'))] # len(df) = 3352
# More filtering criteria? No, by the moment
# Splice variants. I have summarized all variants containing "splice" and checked for their distance to NeasrestSS. df.loc[df['ANNOTATION'].str.contains('splice')].groupby(['ANNOTATION','distNearestSS'])['ID'].count()
# I will filter out variants labeled as "splice_acceptor_variant&splice_region_variant&intron_variant" and "splice_donor_variant&splice_region_variant&intron_variant"
# Variants in those two categories are located further to the SS (as far as -11 and 13 bp)
df['distNearestSS_aux'] = df['distNearestSS']
mask = (df['annonimous_GENE'].str.contains(','))
df_valid = df[mask]
df.loc[mask, 'distNearestSS_aux'] = df_valid.apply(lambda x: x['distNearestSS_aux'].split(',')[int(x['aux'])],axis = 1)

# Filtering splice_acceptor_variant&splice_region_variant&intron_variant in the splicing canonical regions (-2 and -1 bp from the start of the exon)
dfaux = df.loc[df['ANNOTATION'] == 'splice_acceptor_variant&splice_region_variant&intron_variant'] #len(dfaux) = 114
dfaux = dfaux.loc[df['distNearestSS_aux'].astype(float)>-3] #len(dfaux) = 1
df = df.loc[~(df['ANNOTATION'] == 'splice_acceptor_variant&splice_region_variant&intron_variant')] # len(df) = 3238; these variants are further than 3 bp from the canonical splice site
df = pd.concat([df,dfaux]) #len(df) = 3239

# Filteringsplice_donor_variant&splice_region_variant&intron_variant in the splicing canonical regions (+2 and +1 bp from the end of the exon)
dfaux = df.loc[df['ANNOTATION'] == 'splice_donor_variant&splice_region_variant&intron_variant'] #len(dfaux) = 114
dfaux = dfaux.loc[df['distNearestSS_aux'].astype(float)<3] #len(dfaux) = 1
df = df.loc[~(df['ANNOTATION'] == 'splice_donor_variant&splice_region_variant&intron_variant')] # len(df) = 3238; these variants are further than 3 bp from the canonical splice site
df = pd.concat([df,dfaux]) #len(df) = 3127

# Filtering out all variants which are spliceACCEPTOR/DONOR&intron_variant with close distance to NeasrestSS (-3,3)
dfaux = df.loc[~df['ANNOTATION'].str.contains('splice_acceptor_variant&intron_variant')]
dfaux = dfaux.loc[~dfaux['ANNOTATION'].str.contains('splice_donor_variant&intron_variant')]
dff = df.loc[(df['ANNOTATION'].str.contains('splice_acceptor_variant&intron_variant')) & (df['distNearestSS_aux'].astype(float)>-3) ]
dff2 = df.loc[(df['ANNOTATION'].str.contains('splice_donor_variant&intron_variant')) & (df['distNearestSS_aux'].astype(float)<3) ]
df = pd.concat([dfaux,dff, dff2]) #len(df) = 2481

# Create a new column with simplified protein consequence terms
df['ANNOTATION_simple'] = df['ANNOTATION'].apply(lambda x: x.split(' ')[0].split('&')[0].replace('start_lost','start/stop gained/lost').replace('stopgain','start/stop gained/lost').replace('stoploss','start/stop gained/lost').replace('nonframeshift','inframe').split('_')[0])


# Annotate restricted data: per individual NGS location
procdf = pd.read_csv('/path/to/dictionary_proc_pais.csv',sep='|',header = None) # Restricted
procdf_previous = pd.read_csv('/path/to/dictionary_proc.csv',sep='\t') # Restricted
sampledf = pd.read_csv('/path/to/sample_dictionary_20210105.csv',sep='\t',header = None) # Restricted
sample_dict = dict(zip(list(sampledf[1]),list(sampledf[0].astype(str))))
proc_dict = dict(zip(list(procdf[0].astype(str)),list(procdf[1])))
proc_dict_previous = dict(zip(list(procdf_previous['sample'].astype(str)),list(procdf_previous['from'])))
df['annonimous_samples'] = df['samples'].apply(lambda x: ','.join([sample_dict[i] for i in x.split(';')]))
df['annonimous_proc_previous'] = df['annonimous_samples'].apply(lambda x: ','.join([proc_dict_previous[i] for i in x.split(',')]))
df['annonimous_proc'] = df['annonimous_samples'].apply(lambda x: ','.join([proc_dict[i] for i in x.split(',')]))



############## Figure 4A
plt.style.use('default')

df['annonimous_proc_simple'] = df.apply(lambda x: ','.join(list(set(x['annonimous_proc'].split(',')))),axis=1)
df.loc[(df['gnomad_AF'] == 0) & (df['dbSNP_id'] == '.')].groupby('annonimous_proc_simple')['ID'].count()
othercombs =  [x for x in list(set(list(df['annonimous_proc_simple']))) if (x != 'Spain') and (x != 'Colombia') and (x != 'Brazil')]
GENES = ['CACNA1S','CYP2B6','CYP2C9','CYP4F2','DPYD','G6PD','NUDT15','RYR1','SLCO1B1','TPMT','UGT1A1'] # Only ClinAcc PGx genes
df = df.loc[df['SYMBOL'].str.contains('|'.join(GENES))]



f, ax1 = plt.subplots(figsize=(2,3))
dff = df.loc[(df['gnomad_AF'] != 0) | (df['dbSNP_id'] != '.')]
dff_novel = df.loc[(df['gnomad_AF'] == 0) & (df['dbSNP_id'] == '.')]
colordict = dict(zip(['Spain','Colombia','Brazil', '|'.join(othercombs)],['tab:blue','tab:red','tab:green', 'tab:olive']))
bottom = (0,0)
for i in ['Spain','Colombia','Brazil', '|'.join(othercombs)]:
    if i == '|'.join(othercombs):
        labels = 'Mixed'
        countriescount = (len(dff.loc[dff['annonimous_proc_simple'].str.contains(i)]), len(dff_novel.loc[dff_novel['annonimous_proc_simple'].str.contains(i)]))
    else:
        labels = i
        countriescount = (len(dff.loc[dff['annonimous_proc_simple'] == i])),(len(dff_novel.loc[dff_novel['annonimous_proc_simple'] == i]))
    ax1.bar(['Known\nvariants','Novel\nvariants'], countriescount, color = [colordict[i]],label = labels, edgecolor = 'black', bottom = bottom,width = 0.5)#align = 'edge',
    bottom = ( bottom[0] + countriescount[0], bottom[1] + countriescount[1])
ax1.legend(title='Countries', handletextpad=0.5,labelspacing=0.2, fontsize=8,facecolor='white', loc='upper right', bbox_to_anchor=(1.3, 0.8))
plt.ylabel('Count')
plt.style.use('ggplot')
ax1.grid(b=None)#, linestyle='-', linewidth=2)
ax1.set_facecolor('white')
ax1.spines['bottom'].set_color('k')
ax1.spines['top'].set_color('w')
ax1.spines['right'].set_color('w')
ax1.spines['left'].set_color('k')
plt.subplots_adjust(left=0.295, bottom=0.140, right=0.805, top=0.985, wspace=0.17, hspace=0.42)
plt.savefig('/path/to/Figures/Figure_4A.png',format = 'png', dpi = 500)
plt.show()




# Figure 4B (piechart)
f, ax2 = plt.subplots(figsize=(3,2.5))
dfplot = df.loc[(df['gnomad_AF'] == 0) & (df['dbSNP_id'] == '.')].groupby('ANNOTATION_simple')['ID'].count().reset_index();
dfplot.plot(kind='pie',ax=ax2, y='ID', autopct='%1.1f%%',startangle=25, shadow=False, labels = dfplot['ANNOTATION_simple'], labeldistance=np.nan, pctdistance=1.24, textprops={'family':'sans-serif', 'fontsize':9})
ax2.legend(fontsize=8, loc='lower left', bbox_to_anchor=(0.3, 0.1),handletextpad=0.1,labelspacing=0.2)
plt.ylabel('')
nrnovel= len(df.loc[(df['gnomad_AF'] == 0) & (df['dbSNP_id'] == '.')])
ax2.set_title('Novel variants\n' + '(n=' + str(nrnovel) + ')', fontsize=10, pad=-8)
plt.subplots_adjust(right=0.865, bottom=0.0)#, left=0.805, top=0.985, wspace=0.17, hspace=0.42)
plt.savefig('/path/to/Figure_4B.pdf',format = 'pdf', dpi = 500)
plt.show()



#### Figure 4C ####
import numpy as np
df_known = df.loc[~((df['gnomad_AF'] == 0) & (df['dbSNP_id'] == '.'))]
dff_novel = df.loc[(df['gnomad_AF'] == 0) & (df['dbSNP_id'] == '.')]
dff_novel_lof = dff_novel.loc[~dff_novel['ANNOTATION_simple'].str.contains('missense|inframe')]
dff_novel_VUS = dff_novel.loc[dff_novel['ANNOTATION_simple'].str.contains('missense|inframe')]

genesdf = pd.read_csv('../bioMart_transcripts_length.csv',sep='\t')
genesdf = genesdf.loc[genesdf['Gene name'].str.contains('|'.join(GENES))]
genesdf['Transcript length (kb)'] = genesdf['Transcript length (including UTRs and CDS)']/1000

d = dict(dff_novel_lof.loc[dff_novel_lof['SYMBOL'].str.contains('|'.join(GENES))].groupby('SYMBOL')['N_samples'].sum())
genesdf['Nr.Novel LoF variants'] = genesdf['Gene name'].map(d)
genesdf['Nr.Novel Lof vars per kb'] = genesdf['Nr.Novel LoF variants'] / genesdf['Transcript length (kb)']

d = dict(dff_novel_VUS.loc[dff_novel_VUS['SYMBOL'].str.contains('|'.join(GENES))].groupby('SYMBOL')['N_samples'].sum())
genesdf['Nr.Novel VUS variants'] = genesdf['Gene name'].map(d)
genesdf['Nr.Novel VUS vars per kb'] = genesdf['Nr.Novel VUS variants'] / genesdf['Transcript length (kb)']


plotdf = genesdf[['Gene name','Transcript length (kb)','Nr.Novel Lof vars per kb','Nr.Novel LoF variants','Nr.Novel VUS variants','Nr.Novel VUS vars per kb']]

haplo = pd.read_csv('/path/to/haplotypes_20210107.csv',sep='\t')
pheno = pd.read_csv('/path/to/phenotypes_20210107.csv',sep='\t')
genes = list()
values = list()
values_pheno = list()
for gene in GENES:
    if gene == 'CACNA1S':
        genes.append(gene)
        values.append(0)
        values_pheno.append(0)
    elif gene == 'CYP2B6':
        genes.append(gene)
        values.append(len(haplo.loc[~(haplo[gene].isnull())].loc[~haplo[gene].astype(str).str.contains('4|22')]))
        values_pheno.append(len(pheno.loc[~(pheno['Phenotype_' + gene].isnull())]))
    elif gene == 'G6PD':
        genes.append(gene)
        values.append(len(haplo.loc[~(haplo[gene].isnull())]))
        values_pheno.append(len(pheno.loc[pheno['Phenotype_' + gene].astype(str).str.contains('G6PD_Deficient')]))
    else:
        genes.append(gene)
        values.append(len(haplo.loc[~(haplo[gene].isnull())]))
        values_pheno.append(len(pheno.loc[~(pheno['Phenotype_' + gene].isnull())]))
gf = pd.DataFrame()
gf['Gene'] = genes
gf['Nr.ind'] = values_pheno # Considering individuals with actionable phenotypes
gf = gf.sort_values(by='Nr.ind',ascending=False)

d = dict(zip(list(gf['Gene']), list(gf['Nr.ind'])))
plotdf['Nr.actionable alleles'] = plotdf['Gene name'].map(d)
plotdf['PCT Novel Variants'] = [(str(x)+'%').replace('inf%','+36').replace('nan%','0.0%') for x in list((100*plotdf['Nr.Novel LoF variants'] / plotdf['Nr.actionable alleles']).round(decimals=1))]
plotdf.loc[(plotdf['Gene name'].str.contains('CACNA1S|RYR1')),'PCT Novel Variants'] ='NA'
UGT1A1 = np.round(100*(plotdf['Nr.Novel LoF variants'].loc[plotdf['Gene name'] == 'UGT1A1'].values[0] / plotdf['Nr.actionable alleles'].loc[plotdf['Gene name'] == 'UGT1A1'].values[0]), 2)
plotdf.loc[(plotdf['Gene name'].str.contains('UGT1A1')),'PCT Novel Variants'] = str(UGT1A1) + '%'
plt.style.use('ggplot')
plotdf.set_index('Gene name', inplace=True)

f, ax1 = plt.subplots(figsize=(5,3))

plotdf['Nr.Novel LoF variants'].plot(ax=ax1,kind='bar', width=0.75, edgecolor = 'black', label = 'Novel LOF')
# Text on the top of each barplot
labels = list(plotdf['PCT Novel Variants'])
values = list(plotdf['Nr.Novel LoF variants'] + 0.8) #list(plotdf['Nr.Novel LoF variants'])
for i in range(len(plotdf)):
    plt.text(x = i-0.27 , y = values[i]+0.3, s =  labels[i], size = 7)
plt.ylim([0,10])
plt.ylabel('nr. novel LOF variants', fontsize = 10)
ax1.grid(b=None)#, linestyle='-', linewidth=2)
ax1.set_facecolor('white')
ax1.spines['bottom'].set_color('k')
ax1.spines['top'].set_color('w')
ax1.spines['right'].set_color('w')
ax1.spines['left'].set_color('k')
# Adjust the margins
plt.subplots_adjust(bottom= 0.27, top = 0.95, left=0.135, right=0.945)
plt.savefig('/path/to/Figures/Figure_4C.png',format = 'png', dpi = 500)
plt.show()




################################# STATS
# Calculate the mean MAFs of known and novel variants
df_known['MAF'].mean()
dff_novel['MAF'].mean()

# Chi-square test to compare MAFs of known variants versus MAFs of novel variants
from scipy.stats import ttest_ind
ttest_ind(df_known['MAF'],dff_novel['MAF'])

# % Known variants found in only one country
100*(len(df_known.loc[~df_known['annonimous_proc_simple'].str.contains(',')]) / len(df_known))

# % Novel variants found in only one country
100*len(dff_novel.loc[~dff_novel['annonimous_proc_simple'].str.contains(',')]) / len(dff_novel)


# Average % contribution novel LoF variants to actionable alleles:
plotdf['PCT Novel Variants'].str.replace('%','').str.replace('NA','0').astype(float).mean()

# % novel "VUS" (i.e., missense of inframe) in our results
100*len(dff_novel.loc[dff_novel['ANNOTATION_simple'].str.contains('missense|inframe')]) / len(dff_novel)


