import pandas as pd
def multiIDgnomad(CHR, position, ref, alt):
    position_list = str(position).split(',')
    chr_element = ['chr']*len(list(str(position).split(',')))
    chr_list = [str(CHR)]*len(str(position).split(','))
    chr_final = list(map("".join,list(zip(chr_element,chr_list))))
    ref_list = str(ref).split(',')
    alt_list = str(alt).split(',')
    outlist = list(map(":".join, list(zip(chr_final,position_list))))
    outlist = list(map("".join, list(zip(outlist,ref_list))))
    outlist = list(map(">".join, list(zip(outlist,alt_list))))
    outlist = [x.rstrip('>').replace('nan>','').replace('>nan','') for x in outlist]
    outlist = ','.join(outlist)
    return(outlist)

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



df = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
#df['actionable'].loc[(df['SYMBOL'] == 'CYP4F2') & (df['allele'] == '*2')] = 'Yes'
df = df.loc[(df['count_carrier_ids'].astype(str) != 'nan') & (df['actionable'] == 'Yes')].copy()
reference = pd.read_csv('/path/to/reference_HAPLOTYPES_20201130_hg38_hg19.csv',sep='\t')
reference['ID2'] = reference['SYMBOL'] + '_' + reference['allele'].astype(str)
reference['multiID_gnomad'] = reference.apply(lambda x: multiIDgnomad(x['chr'], x['position_hg19_gnomad'],x['ref_gnomad'],x['alt_gnomad']),axis=1)
df['multiID_gnomad'] = df['ID2'].map(dict(zip(list(reference['ID2']),list(reference['multiID_gnomad']))))
sampleinfo = pd.read_csv('/home/jlanillos/CNIO/PGx/dictionary_proc.csv',sep='\t')
sampleinfo['sample'] = sampleinfo['sample'].astype(str)
gender = pd.read_csv('/path/to/gender.csv',sep='\t')
gender_dict = dict(zip(list(gender['annonimousID'].astype(str)),list(gender['gender'])))
country_dict = dict(zip(list(sampleinfo['sample']),list(sampleinfo['from'])))
country_dict_general = dict(zip(list(sampleinfo['sample']),list(sampleinfo['from_general'])))


def calculateNallelesperRS(multiID_gnomad, gts):
    multiID_gnomad = multiID_gnomad.split(',')
    d_al = list()
    for al in multiID_gnomad:
        genotypes = [j.split(',')[multiID_gnomad.index(al)] for j in gts.split(';')]
        heteroz = genotypes.count('0/1')
        homoz = 2*genotypes.count('1/1')
        wt = 2*genotypes.count('0/0') + heteroz # Adding all wt alleles, including those from 0/0 individuals and those from 0/1
        othergenotypes = list()
        for z in list(set(genotypes)):
            if (z != '0/1') and (z != '1/1') and (z != '0/0'):
                othergenotypes.append(z)
                wt = wt + 2*genotypes.count(z) # assume that other values different than 0/1, 1/1 and 0/0 are considered as 0/0 too
        d_al.append(','.join([str(heteroz + homoz),str(wt)]))
    return(';'.join(d_al))

# New function added to solve the count of alleles in X-chrom genes such as G6PD
def correctNallelesperRS_chromX(multiID_gnomad, gts, samples, gender_dict):
    multiID_gnomad = multiID_gnomad.split(',')
    samples = [i.rstrip('_gt') for i in samples.split(',')]
    gender_samples = [gender_dict[i] for i in samples]
    d_al = list()
    for al in multiID_gnomad:
        genotypes = [j.split(',')[multiID_gnomad.index(al)] for j in gts.split(';')]
        genotypes_female = [i for i,j in zip(genotypes, gender_samples) if j == 'F']
        genotypes_male = [i for i,j in zip(genotypes, gender_samples) if j == 'M']
        heteroz = genotypes_female.count('0/1') + genotypes_male.count('0/1')
        homoz = 2*genotypes_female.count('1/1') + genotypes_male.count('1/1')
        wt = 2*genotypes_female.count('0/0') + genotypes_male.count('0/0')
        othergenotypes = list()
        for z in list(set(genotypes)):
            if (z != '0/1') and (z != '1/1') and (z != '0/0'):
                othergenotypes.append(z)
                wt = wt + 2*genotypes_female.count(z) + genotypes_male.count(z) # assume that other values different than 0/1, 1/1 and 0/0 are considered as 0/0 too
        d_al.append(','.join([str(heteroz + homoz),str(wt)]))
    return(';'.join(d_al))

def splitByCountry(carrier_ids, carrier_gt, wt_ids, wt_gt, country_dict, country, key):
    carrier_ids = carrier_ids.replace('_gt','').split(',')
    carrier_gt = carrier_gt.split(';')
    new_carrier_ids = list()
    new_carrier_gt = list()
    wt_ids = wt_ids.replace('_gt','').split(',')
    wt_gt = wt_gt.split(';')
    new_wt_ids = list()
    new_wt_gt = list()
    for i,j in zip(carrier_ids,carrier_gt):
        if country_dict[i] == country:
            new_carrier_ids.append(i+'_gt')
            new_carrier_gt.append(j)
    new_carrier_ids = ','.join(new_carrier_ids)
    new_carrier_gt = ';'.join(new_carrier_gt)
    for i,j in zip(wt_ids,wt_gt):
        if country_dict[i] == country:
            new_wt_ids.append(i+'_gt')
            new_wt_gt.append(j)
    new_wt_ids = ','.join(new_wt_ids)
    new_wt_gt = ';'.join(new_wt_gt)
    if key == 'onlypositivevalues_carrier_ids':
        outlist = new_carrier_ids
    elif key == 'onlypositivevalues_carrier_gt':
        outlist = new_carrier_gt
    elif key == 'onlypositivevalues_wt_ids':
        outlist = new_wt_ids
    elif key == 'onlypositivevalues_wt_gt':
        outlist = new_wt_gt
    return(outlist)


countries = list(set(sampleinfo['from']))
cols = ['onlypositivevalues_carrier_ids','onlypositivevalues_carrier_gt','onlypositivevalues_wt_ids','onlypositivevalues_wt_gt']
for country in countries:
    for col in cols:
        df[country + '_' + col] = df.apply(lambda x: splitByCountry(x['onlypositivevalues_carrier_ids'],x['onlypositivevalues_carrier_gt'],x['onlypositivevalues_wt_ids'],x['onlypositivevalues_wt_gt'], country_dict, country, col),axis = 1)
    df['aux'] = df[country + '_' +'onlypositivevalues_carrier_gt'].astype(str) + ';' + df[country + '_' + 'onlypositivevalues_wt_gt'].astype(str)
    df['aux'] = df['aux'].apply(lambda x: x.replace('nan;','').replace(';nan','').rstrip(';').lstrip(';'))
    df['aux_ids'] = df[country + '_' +'onlypositivevalues_carrier_ids'].astype(str) + ',' + df[country + '_' + 'onlypositivevalues_wt_ids'].astype(str)
    df['aux_ids'] = df['aux_ids'].apply(lambda x: x.replace('nan,','').replace(',nan','').rstrip(',').lstrip(','))
    df[country + '_N_alleles'] = df.apply(lambda x: calculateNallelesperRS(x['multiID_gnomad'], x['aux']) if str(x['chr']) != 'X' else correctNallelesperRS_chromX(x['multiID_gnomad'], x['aux'], x['aux_ids'], gender_dict),axis=1)

countries = ['LATAM'] # list(set(sampleinfo['from_general']))
cols = ['onlypositivevalues_carrier_ids','onlypositivevalues_carrier_gt','onlypositivevalues_wt_ids','onlypositivevalues_wt_gt']
for country in countries:
    for col in cols:
        df[country + '_' + col] = df.apply(lambda x: splitByCountry(x['onlypositivevalues_carrier_ids'],x['onlypositivevalues_carrier_gt'],x['onlypositivevalues_wt_ids'],x['onlypositivevalues_wt_gt'], country_dict_general, country, col),axis = 1)
    df['aux'] = df[country + '_' +'onlypositivevalues_carrier_gt'].astype(str) + ';' + df[country + '_' + 'onlypositivevalues_wt_gt'].astype(str)
    df['aux'] = df['aux'].apply(lambda x: x.replace('nan;','').replace(';nan','').rstrip(';').lstrip(';'))
    df['aux_ids'] = df[country + '_' +'onlypositivevalues_carrier_ids'].astype(str) + ',' + df[country + '_' + 'onlypositivevalues_wt_ids'].astype(str)
    df['aux_ids'] = df['aux_ids'].apply(lambda x: x.replace('nan,','').replace(',nan','').rstrip(',').lstrip(','))
    df[country + '_N_alleles'] = df.apply(lambda x: calculateNallelesperRS(x['multiID_gnomad'], x['aux']) if str(x['chr']) != 'X' else correctNallelesperRS_chromX(x['multiID_gnomad'], x['aux'], x['aux_ids'], gender_dict),axis=1)

country = 'ALL'
df['aux'] = df['onlypositivevalues_carrier_gt'].astype(str) + ';' + df['onlypositivevalues_wt_gt'].astype(str)
df['aux'] = df['aux'].apply(lambda x: x.replace('nan;','').replace(';nan','').rstrip(';').lstrip(';'))
df['aux_ids'] = df['onlypositivevalues_carrier_ids'].astype(str) + ',' + df['onlypositivevalues_wt_ids'].astype(str)
df['aux_ids'] = df['aux_ids'].apply(lambda x: x.replace('nan,','').replace(',nan','').rstrip(',').lstrip(','))
df[country + '_N_alleles'] = df.apply(lambda x: calculateNallelesperRS(x['multiID_gnomad'], x['aux']) if str(x['chr']) != 'X' else correctNallelesperRS_chromX(x['multiID_gnomad'], x['aux'], x['aux_ids'], gender_dict),axis=1)


gnomad = pd.read_csv('/path/to/gnomAD_info/gnomAD_rs_alleles_pop.info.csv',sep='\t') # Done by 20210117_getVariantInfo2jsonandparsetoref.py
cols = ['_'.join(x.split('_')[1:]) for x in list(gnomad.columns) if 'an_' in x]
cols = [''] + cols # To calculate gnomad total frequency too
for col in cols:
    if col == '':
        gnomad['N_alleles' +'_' + col] = gnomad['an'] *gnomad['af']
        gnomad['N_alleles' +'_' + col] = gnomad['N_alleles' + '_' + col].round()
        gnomad['N_alleles_wt' +'_' + col] = gnomad['an'] - gnomad['N_alleles' + '_' + col]
    else:
        gnomad['N_alleles' +'_' + col] = gnomad['an_' + col] *gnomad['af_' + col]
        gnomad['N_alleles' +'_' + col] = gnomad['N_alleles' + '_' + col].round()
        gnomad['N_alleles_wt' +'_' + col] = gnomad['an_' + col] - gnomad['N_alleles' + '_' + col]


ref_alleles = pd.read_csv('/path/to/gnomAD_info/allele2VarId.csv',sep='\t') # This s list of alleles with uniques rs which cover all alleles with results, curated to be used as template in the next step: extract the N_alleles count from each country
# for those chosen alleles/rs(s). We have to calculate MERGED_GVCFs.csv again and modify next step to avoid this "manual" need
ref_alleles['aux'] = ref_alleles['allele']  + '_' + ref_alleles['id']
for country in list(set(sampleinfo['from'])) + ['LATAM', 'ALL']:
    p = ';'.join(list(df[country + '_N_alleles'])).split(';')
    m = ','.join(list(df['multiID_gnomad'])).split(',')
    al = list()
    for i,j in zip(list(df['ID2']), list(df['multiID_gnomad'].apply(lambda x: len(x.split(','))))):
        al = al + [i]*j
    gf = pd.DataFrame({'allele': al, 'id': m, 'count': p})
    gf['allele'] = gf['allele'] + '_' + gf['id']
    d = dict(zip(list(gf['allele']),list(gf['count'].apply(lambda x: int(x.split(',')[0])))))
    ref_alleles['N_alleles_' + country] = ref_alleles['aux'].map(d)
    d = dict(zip(list(gf['allele']),list(gf['count'].apply(lambda x: int(x.split(',')[1])))))
    ref_alleles['N_alleles_wt_' + country] = ref_alleles['aux'].map(d)

cols = [x for x in list(ref_alleles.columns) if 'N_alleles_' in x]
for col in cols:
    d = dict(zip(list(ref_alleles['id']),list(ref_alleles[col])))
    gnomad[col] = gnomad['fields'].map(d)

##gnomad.to_csv('gnomAD_info/ResultsAnnot_gnomAD_rs_alleles_pop.info.csv',sep='\t')
gnomad.to_csv('path/to/gnomAD_info/ResultsAnnot_gnomAD_rs_alleles_pop.info_splittedbycountry.csv',sep='\t')


import numpy as np
from scipy.stats import chi2_contingency
##gnomad = pd.read_csv('gnomAD_info/ResultsAnnot_gnomAD_rs_alleles_pop.info.csv',sep='\t')
gnomad = pd.read_csv('path/to/gnomAD_info/ResultsAnnot_gnomAD_rs_alleles_pop.info_splittedbycountry.csv',sep='\t')

pvalthre = 0.05
ourdatagroups = ['Spain', 'Colombia','Brazil', 'LATAM', 'ALL']
gnomadcomp = ['nfe', 'amr']
gnomadchi2df = gnomad.copy() #gnomad[['SYMBOL','alleles','fields','rsid']]

for ourdata in ourdatagroups:
    for gnomadgroup in gnomadcomp:
        chi2 = gnomad.loc[gnomad['N_alleles_' + ourdata]>1].apply(lambda x: chi2_contingency(np.array([[x['N_alleles_' + ourdata],x['N_alleles_wt_' + ourdata]],[x['N_alleles_' + gnomadgroup],x['N_alleles_wt_' + gnomadgroup]]]))[1],axis = 1).values
        fields =  gnomad['fields'].loc[gnomad['N_alleles_' + ourdata]>1].values
        correction = pvalthre / len(chi2)
        dff = pd.DataFrame({'fields':fields , 'chi2':chi2})
        dff = dff.loc[dff['chi2'] <= correction]
        d = dict(zip(list(dff['fields']),list(dff['chi2'])))
        gnomadchi2df['chi2pval_' + ourdata + '_' + gnomadgroup] = gnomadchi2df['fields'].map(d)
    # Compare to gnomad_AF
    chi2 = gnomad.loc[gnomad['N_alleles_' + ourdata]>1].apply(lambda x: chi2_contingency(np.array([[x['N_alleles_' + ourdata],x['N_alleles_wt_' + ourdata]],[x['N_alleles_'],x['N_alleles_wt_']]]))[1],axis = 1).values
    fields =  gnomad['fields'].loc[gnomad['N_alleles_' + ourdata]>1].values
    correction = pvalthre / len(chi2)
    dff = pd.DataFrame({'fields':fields , 'chi2':chi2})
    dff = dff.loc[dff['chi2'] <= correction]
    d = dict(zip(list(dff['fields']),list(dff['chi2'])))
    gnomadchi2df['chi2pval_' + ourdata + '_' + 'GENERAL'] = gnomadchi2df['fields'].map(d)


#gnomadchi2df.loc[~gnomadchi2df['chi2pval_ESPAÑA_nfe'].isnull()].sort_values(by=['SYMBOL'])
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
fields = list()
#fields = list(gnomadchi2df['fields'].loc[(~gnomadchi2df['chi2pval_ALL_GENERAL'].isnull())])
fields_Spain = list(gnomadchi2df['fields'].loc[(~gnomadchi2df['chi2pval_Spain_nfe'].isnull())])
fields_LATAM = list(gnomadchi2df['fields'].loc[(~gnomadchi2df['chi2pval_LATAM_amr'].isnull())])
fields_Colombia = list(gnomadchi2df['fields'].loc[(~gnomadchi2df['chi2pval_Colombia_amr'].isnull())])
fields_Brazil = list(gnomadchi2df['fields'].loc[(~gnomadchi2df['chi2pval_Brazil_amr'].isnull())])

fields = fields +  fields_Spain + fields_LATAM + fields_Colombia + fields_Brazil
fields = list(set(fields))

chi2plot = gnomadchi2df.loc[gnomadchi2df['fields'].str.contains('|'.join(fields))]
chi2plot = chi2plot.sort_values(['SYMBOL', 'fields'], ascending=[True, True])
chi2plot['af_Spain'] = chi2plot['N_alleles_Spain'] / (chi2plot['N_alleles_Spain']  + chi2plot['N_alleles_wt_Spain'])
chi2plot['af_ALL'] = chi2plot['N_alleles_ALL'] / (chi2plot['N_alleles_ALL']  + chi2plot['N_alleles_wt_ALL'])
chi2plot['af_LATAM'] = chi2plot['N_alleles_LATAM'] / (chi2plot['N_alleles_LATAM']  + chi2plot['N_alleles_wt_LATAM'])
chi2plot['af_Colombia'] = chi2plot['N_alleles_Colombia'] / (chi2plot['N_alleles_Colombia']  + chi2plot['N_alleles_wt_Colombia'])
chi2plot['af_Brazil'] = chi2plot['N_alleles_Brazil'] / (chi2plot['N_alleles_Brazil']  + chi2plot['N_alleles_wt_Brazil'])

chi2plot['full_fields'] = chi2plot['SYMBOL'] + ' - ' + chi2plot['rsid']# + '\n' + '(' + chi2plot['fields'] + ')'
chi2AUX = chi2plot.copy()


# Keep only chi2 comparison which make sense: Spain vs NFE, Colombia vs AMR, Brazil vs AMR, LATAM vs AMR and All vs gnomad_general
cols2remove = ['chi2pval_Spain_amr', 'chi2pval_Spain_GENERAL', 'chi2pval_Colombia_nfe', 'chi2pval_Colombia_GENERAL', 'chi2pval_Brazil_nfe','chi2pval_Brazil_GENERAL', 'chi2pval_LATAM_nfe', 'chi2pval_LATAM_GENERAL', 'chi2pval_ALL_nfe', 'chi2pval_ALL_amr']
chi2plotaux = chi2plot.drop(cols2remove, axis = 1)
chi2plot = chi2plotaux.copy()



f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10,5), gridspec_kw={'width_ratios':[2,3,5]})
freqs_pops = ['af','af_amr','af_nfe','af_Spain','af_ALL','af_LATAM','af_Colombia','af_Brazil']

chi2plot = chi2plotaux.loc[chi2plotaux[freqs_pops].mean(axis=1) >= 0.1]
ax1.scatter(chi2plot['full_fields'], chi2plot['af'], c = 'black', s = 30, label = 'gnomAD_af', marker = 's', edgecolor ="black")#, transform=trans+offset(2))
ax1.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(list(set(fields_LATAM + fields_Brazil + fields_Colombia))))], chi2plot['af_amr'].loc[chi2plot['fields'].str.contains('|'.join(list(set(fields_LATAM + fields_Brazil + fields_Colombia))))],alpha=0.8, c = 'blue', s = 30, label = 'gnomad_amr', marker = 's', edgecolor ="blue")
#ax1.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))], chi2plot['af_amr'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))],alpha=0.8, c = 'blue', s = 30, label = 'gnomad_amr', marker = 's', edgecolor ="blue")
ax1.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], chi2plot['af_nfe'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], c = 'red',alpha=0.8, s = 30, label = 'gnomad_nfe', marker = 's', edgecolor ="red")
ax1.scatter(chi2plot['full_fields'], chi2plot['af_ALL'], c = 'white', s = 30, label = 'af_ALL', marker = 'o', edgecolor ="black")
ax1.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))], chi2plot['af_LATAM'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))],alpha=0.8, c = 'white', s = 30, label = 'af_LATAM', marker = 'o', edgecolor ="blue")
ax1.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], chi2plot['af_Spain'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], c = 'white',alpha=0.8, s = 30, label = 'af_Spain', marker = 'o', edgecolor ="red")
ax1.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Colombia))], chi2plot['af_Colombia'].loc[chi2plot['fields'].str.contains('|'.join(fields_Colombia))],alpha=0.8, c = 'white', s = 30, label = 'af_Colombia', marker = 'o', edgecolor ="green")
ax1.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Brazil))], chi2plot['af_Brazil'].loc[chi2plot['fields'].str.contains('|'.join(fields_Brazil))],alpha=0.8, c = 'white', s = 30, label = 'af_Brazil', marker = 'o', edgecolor ="purple")
#ax1.xticks(fontsize = 8,rotation=90)
ax1.tick_params(axis='x',size = 8,labelrotation=90)
ax1.set_ylabel('Frequency',size=12)

chi2plot = chi2plotaux.loc[(chi2plotaux[freqs_pops].mean(axis=1) < 0.1) & (chi2plotaux[freqs_pops].mean(axis=1) >= 0.01)]
ax2.scatter(chi2plot['full_fields'], chi2plot['af'], c = 'black', s = 30, label = 'gnomAD_af', marker = 's', edgecolor ="black")#, transform=trans+offset(2))
ax2.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(list(set(fields_Brazil + fields_LATAM + fields_Colombia))))], chi2plot['af_amr'].loc[chi2plot['fields'].str.contains('|'.join(list(set(fields_Brazil + fields_LATAM + fields_Colombia))))],alpha=0.8, c = 'blue', s = 30, label = 'gnomad_amr', marker = 's', edgecolor ="blue")
ax2.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], chi2plot['af_nfe'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], c = 'red',alpha=0.8, s = 30, label = 'gnomad_nfe', marker = 's', edgecolor ="red")
ax2.scatter(chi2plot['full_fields'], chi2plot['af_ALL'], c = 'white', s = 30, label = 'af_ALL', marker = 'o', edgecolor ="black")
ax2.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))], chi2plot['af_LATAM'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))],alpha=0.8, c = 'white', s = 30, label = 'af_LATAM', marker = 'o', edgecolor ="blue")
ax2.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], chi2plot['af_Spain'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], c = 'white',alpha=0.8, s = 30, label = 'af_Spain', marker = 'o', edgecolor ="red")
ax2.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Colombia))], chi2plot['af_Colombia'].loc[chi2plot['fields'].str.contains('|'.join(fields_Colombia))],alpha=0.8, c = 'white', s = 30, label = 'af_Colombia', marker = 'o', edgecolor ="green")
ax2.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Brazil))], chi2plot['af_Brazil'].loc[chi2plot['fields'].str.contains('|'.join(fields_Brazil))],alpha=0.8, c = 'white', s = 30, label = 'af_Brazil', marker = 'o', edgecolor ="purple")
ax2.tick_params(axis='x',size = 8,labelrotation=90)
#ax2.xticks(fontsize = 8,rotation=90)

chi2plot = chi2plotaux.loc[(chi2plotaux[freqs_pops].mean(axis=1) < 0.01)]
ax3.scatter(chi2plot['full_fields'], chi2plot['af'], c = 'black', s = 30, label = 'gnomAD_af', marker = 's', edgecolor ="black")#, transform=trans+offset(2))
ax3.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))], chi2plot['af_amr'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))],alpha=0.8, c = 'blue', s = 30, label = 'gnomad_amr', marker = 's', edgecolor ="blue")
ax3.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], chi2plot['af_nfe'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], c = 'red',alpha=0.8, s = 30, label = 'gnomad_nfe', marker = 's', edgecolor ="red")
ax3.scatter(chi2plot['full_fields'], chi2plot['af_ALL'], c = 'white', s = 30, label = 'af_ALL', marker = 'o', edgecolor ="black")
ax3.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))], chi2plot['af_LATAM'].loc[chi2plot['fields'].str.contains('|'.join(fields_LATAM))],alpha=0.8, c = 'white', s = 30, label = 'af_LATAM', marker = 'o', edgecolor ="blue")
ax3.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], chi2plot['af_Spain'].loc[chi2plot['fields'].str.contains('|'.join(fields_Spain))], c = 'white',alpha=0.8, s = 30, label = 'af_Spain', marker = 'o', edgecolor ="red")
ax3.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Colombia))], chi2plot['af_Colombia'].loc[chi2plot['fields'].str.contains('|'.join(fields_Colombia))],alpha=0.8, c = 'white', s = 30, label = 'af_Colombia', marker = 'o', edgecolor ="green")
ax3.scatter(chi2plot['full_fields'].loc[chi2plot['fields'].str.contains('|'.join(fields_Brazil))], chi2plot['af_Brazil'].loc[chi2plot['fields'].str.contains('|'.join(fields_Brazil))],alpha=0.8, c = 'white', s = 30, label = 'af_Brazil', marker = 'o', edgecolor ="purple")
ax3.tick_params(axis='x',size = 8,labelrotation=90)
#ax2.xticks(fontsize = 8,rotation=90)




plt.legend(handletextpad=0.05,labelspacing=0.1, bbox_to_anchor=(0.6,0.5))
plt.subplots_adjust(left=0.060, bottom=0.365, right=0.995, top=0.970, wspace=0.20, hspace=0.20)
#plt.show()
# pval: 0.00139
##f.suptitle('Spanish vs gnomAD NFE or LatAm countries vs gnomAD AMR (Chi-square p-val = ' + str(np.round(correction, decimals=4)) + ')', fontsize=11)
#plt.savefig('path/to/Figures/Figure_2D_splittedbyLatamCountry.png',format = 'png', dpi = 1000)
plt.savefig('0_PGx_manuscript/CPT/Figures/Figure_2D_splittedbyLatamCountry.png',format = 'png', dpi = 1000)
