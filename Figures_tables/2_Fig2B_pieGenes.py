# Script to plot Figure 2B

df = pd.read_csv('/path/to/haplotypes_20210107.csv',sep='\t')
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
alleles['actionable'].loc[(alleles['SYMBOL'] == 'CYP4F2') & (alleles['allele'] == '*2')] = 'Yes'
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))
genes = list()
values = list()
for gene in GENES:
    genes.append(gene)
    values.append(len(df.loc[~(df[gene].isnull())]))
gf = pd.DataFrame()
gf['Gene'] = genes
gf['Nr.ind'] = values
gf = gf.sort_values(by='Nr.ind',ascending=False)

f, ax1 = plt.subplots(figsize=(6,5))
gf.plot(kind='pie',ax=ax1, y='Nr.ind', autopct='%1.1f%%',startangle=0, explode = [0.01, 0.01, 0.01, 0.01, 0.01, 0.03, 0.06, 0.1, 0.3, 0.60 ], shadow=False, labels = gf['Gene'], labeldistance=np.nan, pctdistance=1.20, textprops={'family':'sans-serif', 'fontsize':16})
ax1.legend(fontsize=12, loc='upper left',handletextpad=0.15,labelspacing=0.2, markerscale=0.2,bbox_to_anchor=(0.93, 1.05))
plt.ylabel('')
plt.tight_layout()
plt.subplots_adjust(right=0.782) #top=0.95, wspace=0.14, hspace=0.24,left=0.15, bottom=0.080
plt.savefig('/path/to/Figures/Figure_2B_pie_nrindpergene.png',format = 'png', dpi = 500)
plt.show()
