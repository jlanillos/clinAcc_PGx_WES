# Script uised to plot Figure 2C
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
df = pd.read_csv('/path/to/haplotypes_20210107.csv',sep='\t')
alleles = pd.read_csv('/path/to/Alleles_20201228.csv',sep='\t')
alleles = alleles.loc[(alleles['count_carrier_ids'].astype(str) != 'nan') & (alleles['actionable'] == 'Yes')].copy()
GENES = list(set(list(alleles['SYMBOL'])))
GENES.sort()

ncols = 2
nrows = int(np.ceil( len(GENES) / ncols))
x3 = list()
for i in list(np.arange(0, nrows, 1)):
    x3.append((0,i))
    x3.append((1,i))
d_coords = dict(zip(GENES, x3))

def countdiffalleles(x):
    count = 0
    for i in x.split(','):
        i = i.split('_')[1]
        if '*1/' in i:
            i = i.split('/')
            i.remove('*1')
        else:
            i = i.split('/')
        count = count + len(i)
    return(count)



def prepareDF4plot(df, country, GENE, RF):
    #d = dict()
    d_comphet = dict()
    otherallelelist = list()
    dff = df.loc[df['from'].str.contains(country)].groupby(GENE)['sample'].count().reset_index()
    dff = dff.loc[~(dff[GENE] == '')].copy()
    nr_found = (dff[GENE].apply(lambda x: countdiffalleles(x)) * dff['sample']).sum()
    posalleles = list(set([i.split('_')[1].split(',')[0].split('/')[1] for i in list(set(','.join(list(dff[GENE])).split(',')))]))
    allele_names = list()
    allele_counts = list()
    allele_counts_aux = list()
    allele_names_comphet = list()
    allele_counts_comphet = list()
    allele_counts_aux_comphet = list()
    for al in posalleles:
        posalleles_aux = [x for x in posalleles if x != al]
        homozall = [i for i in list(dff[GENE]) if '/'.join([al,al]) in i]
        heterozall = [i for i in list(dff[GENE]) if (al in i) and not ('/'.join([al,al]) in i)] #[i for i in list(contigency.columns) if '/'.join(['*1',al]) in i]
        comphet = [i for i in list(dff[GENE]) if (al in i) and not ('/'.join([al,al]) in i) and  any(x in i for x in posalleles_aux)]
        dff['aux'] = dff[GENE].str.replace('*','.')
        if homozall == []:
            nr_homoz = 0
        else:
            nr_homoz = dff['sample'].loc[dff['aux'].str.contains('|'.join([i.replace('*','.') for i in homozall]))].sum()
        if heterozall == []:
            nr_heteroz = 0
        else:
            nr_heteroz = dff['sample'].loc[dff['aux'].str.contains('|'.join([i.replace('*','.') for i in heterozall]))].sum()
        if comphet == []:
            nr_comhet = 0
        else:
            nr_comhet = dff['sample'].loc[dff['aux'].str.contains('|'.join([i.replace('*','.') for i in comphet]))].sum() #nr_heteroz -
        count = 100*(nr_heteroz + nr_homoz*2) / nr_found
        count_comphet = 100*(nr_comhet + nr_homoz*2) / nr_found
        count_nocomphet = count - count_comphet
        count_aux =  (nr_heteroz + nr_homoz*2)
        #d[al] = [count]
        allele_names.append(al)
        allele_counts.append(count)
        allele_names_comphet.append(al + '_COMPHET')
        allele_names_comphet.append(al + '_NOCOMPHET')
        allele_counts_comphet.append(count_comphet)
        allele_counts_comphet.append(count_nocomphet)
        allele_counts_aux.append(count_aux)
    rf = pd.DataFrame({'alleles':allele_names,'count':allele_counts})
    rf_comphet = pd.DataFrame({'alleles':allele_names_comphet,'count':allele_counts_comphet})
    if GENE == 'DPYD':
        if country == 'Ecuador|EUROPA|-|Chile|Brazil|USA|Mexico|Portugal|Argentina|Spain|Colombia':
            n_alleles_aux = 5
        else:
            n_alleles_aux = 4
        d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[0:n_alleles_aux]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[0:n_alleles_aux]['count'])]))
        for w in list(d.keys()):
            d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
            d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]

            otheralleleslength = len(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:])
            otherallelelist = list(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['alleles'].values)
            if otheralleleslength > 1:
                #d['Other alleles ' + '(N=' + str(otheralleleslength) + ')'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                d['Other alleles'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['count'].sum()]
                d_comphet['Other alleles_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet['Other alleles_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
            elif otheralleleslength == 1:
                theother = rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['alleles'].values[0]
                d[theother] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['count'].sum()]
                d_comphet[theother +'_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet[theother +'_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
        genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
        genecountdf_comhet = pd.DataFrame(d_comphet)
    elif GENE == 'CYP2B6':
        n_alleles_aux = 3
        d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[0:n_alleles_aux]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[0:n_alleles_aux]['count'])]))
        for w in list(d.keys()):
            d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
            d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]

            otheralleleslength = len(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:])
            otherallelelist = list(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['alleles'].values)
            if otheralleleslength > 1:
                #d['Other alleles ' + '(N=' + str(otheralleleslength) + ')'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                d['Other alleles'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['count'].sum()]
                d_comphet['Other alleles_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet['Other alleles_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
            elif otheralleleslength == 1:
                theother = rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['alleles'].values[0]
                d[theother] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['count'].sum()]
                d_comphet[theother +'_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet[theother +'_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
        genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
        genecountdf_comhet = pd.DataFrame(d_comphet)
    elif GENE == 'CYP2C9':
        n_alleles_aux = 2
        d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[0:n_alleles_aux]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[0:n_alleles_aux]['count'])]))
        for w in list(d.keys()):
            d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
            d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]

            otheralleleslength = len(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:])
            otherallelelist = list(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['alleles'].values)
            if otheralleleslength > 1:
                #d['Other alleles ' + '(N=' + str(otheralleleslength) + ')'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                d['Other alleles'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['count'].sum()]
                d_comphet['Other alleles_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet['Other alleles_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
            elif otheralleleslength == 1:
                theother = rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['alleles'].values[0]
                d[theother] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['count'].sum()]
                d_comphet[theother +'_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet[theother +'_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
        genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
        genecountdf_comhet = pd.DataFrame(d_comphet)
    elif GENE == 'UGT1A1':
        d = dict(zip(list(rf.sort_values(by=['count'], ascending=False)['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False)['count'])]))
        for w in list(d.keys()):
            d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
            d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]
        genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
        genecountdf_comhet = pd.DataFrame(d_comphet)
    elif GENE == 'RYR1':
        d = dict(zip(list(rf.sort_values(by=['count'], ascending=False)['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False)['count'])]))
        for w in list(d.keys()):
            d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
            d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]
        genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
        genecountdf_comhet = pd.DataFrame(d_comphet)
    elif GENE == 'TPMT':
        d = dict(zip(list(rf.sort_values(by=['count'], ascending=False)['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False)['count'])]))
        for w in list(d.keys()):
            d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
            d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]
        genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
        genecountdf_comhet = pd.DataFrame(d_comphet)
    else:
        if len(posalleles) > 2:
            d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[[0,1]]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[[0,1]]['count'])]))
            for w in list(d.keys()):
                d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
                d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]

                otheralleleslength = len(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:])
                otherallelelist = list(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['alleles'].values)
                if otheralleleslength > 1:
                    #d['Other alleles ' + '(N=' + str(otheralleleslength) + ')'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                    d['Other alleles'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                    d_comphet['Other alleles_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                    d_comphet['Other alleles_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
                else:
                    theother = rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['alleles'].values[0]
                    d[theother] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                    d_comphet[theother +'_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                    d_comphet[theother +'_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
            genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
            genecountdf_comhet = pd.DataFrame(d_comphet)
                    #genecountdf['from'] = country
        elif len(posalleles) == 2:
            d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[[0,1]]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[[0,1]]['count'])]))
            genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
            for w in list(d.keys()):
                d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','.')) == w.replace('*','.')+'_COMPHET'].sum()]
                d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','.')) == w.replace('*','.')+'_NOCOMPHET'].sum()]
            genecountdf_comhet = pd.DataFrame(d_comphet)
        #genecountdf['from'] = country
        elif len(posalleles) == 1:
            d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[[0]]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[[0]]['count'])]))
            for w in list(d.keys()):
                d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','.')) == w.replace('*','.')+'_COMPHET'].sum()]
                d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','.')) == w.replace('*','.')+'_NOCOMPHET'].sum()]
            genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
            genecountdf_comhet = pd.DataFrame(d_comphet)
            #genecountdf['from'] = country

    rf_aux = pd.DataFrame({'Gene':[GENE]*len(allele_names),'alleles':allele_names,'count':allele_counts_aux, 'country': [country]*len(allele_counts_aux), 'nr_found': [nr_found]*len(allele_counts_aux)})
    return genecountdf, genecountdf_comhet, otherallelelist, rf_aux



def prepareDF4plot_chrX(df, country, GENE, RF):
    #d = dict()
    d_comphet = dict()
    otherallelelist = list()
    dff = df.loc[df['from'].str.contains(country)].groupby(GENE)['sample'].count().reset_index()
    dff = dff.loc[~(dff[GENE] == '')].copy()
    dff_females = df.loc[(df['from'].str.contains(country)) & (df['gender'] == 'F')].groupby(GENE)['sample'].count().reset_index()
    dff_females = dff_females.loc[~(dff_females[GENE] == '')].copy()
    dff_males = df.loc[(df['from'].str.contains(country)) & (df['gender'] == 'M')].groupby(GENE)['sample'].count().reset_index()
    dff_males = dff_males.loc[~(dff_males[GENE] == '')].copy()
    nr_found_females = (dff_females[GENE].apply(lambda x: countdiffalleles(x)) * dff_females['sample']).sum()
    nr_found_males = (dff_males[GENE].apply(lambda x: countdiffalleles(x)) * dff_males['sample'] / 2).sum() # or dff_males['sample'].sum()
    nr_found = nr_found_females + nr_found_males
    posalleles = list(set([i.split('_')[1].split(',')[0].split('/')[1] for i in list(set(','.join(list(dff[GENE])).split(',')))]))
    allele_names = list()
    allele_counts = list()
    allele_counts_aux = list()
    allele_names_comphet = list()
    allele_counts_comphet = list()
    allele_counts_aux_comphet = list()
    for al in posalleles:
        posalleles_aux = [x for x in posalleles if x != al]
        homozall = [i for i in list(dff[GENE]) if '/'.join([al,al]) in i]
        heterozall = [i for i in list(dff[GENE]) if (al in i) and not ('/'.join([al,al]) in i)] #[i for i in list(contigency.columns) if '/'.join(['*1',al]) in i]
        comphet = [i for i in list(dff[GENE]) if (al in i) and not ('/'.join([al,al]) in i) and  any(x in i for x in posalleles_aux)]
        dff['aux'] = dff[GENE].str.replace('*','.')
        dff_females['aux'] = dff_females[GENE].str.replace('*','.')
        dff_males['aux'] = dff_males[GENE].str.replace('*','.')
        if homozall == []:
            #nr_homoz = 0
            nr_homoz_female = 0
            nr_homoz_male = 0
        else:
            #nr_homoz = dff['sample'].loc[dff['aux'].str.contains('|'.join([i.replace('*','.') for i in homozall]))].sum()
            nr_homoz_female = dff_females['sample'].loc[dff_females['aux'].str.contains('|'.join([i.replace('*','.') for i in homozall]))].sum()
            nr_homoz_male = dff_males['sample'].loc[dff_males['aux'].str.contains('|'.join([i.replace('*','.') for i in homozall]))].sum()
        if heterozall == []:
            nr_heteroz = 0
        else:
            nr_heteroz = dff['sample'].loc[dff['aux'].str.contains('|'.join([i.replace('*','.') for i in heterozall]))].sum() # Assumption: there are not heterozygous males (checked genotypes before)
        if comphet == []:
            nr_comhet = 0
        else:
            nr_comhet = dff['sample'].loc[dff['aux'].str.contains('|'.join([i.replace('*','.') for i in comphet]))].sum() #nr_heteroz -
        nr_homoz = 2*nr_homoz_female + nr_homoz_male # Number of alleles coming from homozygous individuals, taking into account that male individuals are hemizygous, an thus do add one allele
        count = 100*(nr_heteroz + nr_homoz) / nr_found
        count_comphet = 100*(nr_comhet + nr_homoz) / nr_found
        count_nocomphet = count - count_comphet
        count_aux =  (nr_heteroz + nr_homoz)
        #d[al] = [count]
        allele_names.append(al)
        allele_counts.append(count)
        allele_names_comphet.append(al + '_COMPHET')
        allele_names_comphet.append(al + '_NOCOMPHET')
        allele_counts_comphet.append(count_comphet)
        allele_counts_comphet.append(count_nocomphet)
        allele_counts_aux.append(count_aux)
    rf = pd.DataFrame({'alleles':allele_names,'count':allele_counts})
    rf_comphet = pd.DataFrame({'alleles':allele_names_comphet,'count':allele_counts_comphet})
    n_alleles_aux = 7
    if GENE == 'G6PD':
        d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[0:n_alleles_aux]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[0:n_alleles_aux]['count'])]))
        for w in list(d.keys()):
            d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
            d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]

            otheralleleslength = len(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:])
            otherallelelist = list(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['alleles'].values)
            if otheralleleslength > 1:
                #d['Other alleles ' + '(N=' + str(otheralleleslength) + ')'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                d['Other alleles'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['count'].sum()]
                d_comphet['Other alleles_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet['Other alleles_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
            elif otheralleleslength == 1:
                theother = rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['alleles'].values[0]
                d[theother] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[n_alleles_aux:]['count'].sum()]
                d_comphet[theother +'_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet[theother +'_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
        genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
        genecountdf_comhet = pd.DataFrame(d_comphet)
        genecountdf = genecountdf.rename(columns={'Seattle Lodi Modena FerraraII Athens-like':'Seattle*'})
        genecountdf = genecountdf.rename(columns={'Mediterranean Dallas Panama‚Sassari Cagliari Birmingham':'Mediterranean*'})
        genecountdf = genecountdf.rename(columns={'G6PDA-968C-376G':'A-968C-376G'})
        genecountdf = genecountdf.rename(columns={'Union Maewo Chinese-2 Kalo':'Union*'})
        genecountdf_comhet = genecountdf_comhet.rename(columns={'Seattle Lodi Modena FerraraII Athens-like_COMPHET':'Seattle*_COMPHET','Seattle Lodi Modena FerraraII Athens-like_NOCOMPHET':'Seattle*_NOCOMPHET'})
        genecountdf_comhet = genecountdf_comhet.rename(columns={'Mediterranean Dallas Panama‚Sassari Cagliari Birmingham_COMPHET':'Mediterranean*_COMPHET','Mediterranean Dallas Panama‚Sassari Cagliari Birmingham_NOCOMPHET':'Mediterranean*_NOCOMPHET'})
        genecountdf_comhet = genecountdf_comhet.rename(columns={'G6PDA-968C-376G_COMPHET':'A-968C-376G_COMPHET','G6PDA-968C-376G_NOCOMPHET':'A-968C-376G_NOCOMPHET'})
        genecountdf_comhet = genecountdf_comhet.rename(columns={'Union Maewo Chinese-2 Kalo_COMPHET':'Union*_COMPHET','Union Maewo Chinese-2 Kalo_NOCOMPHET':'Union*_NOCOMPHET'})
    else:
        if len(posalleles) > n_alleles_aux:
            d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[[0,1]]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[[0,1]]['count'])]))
            for w in list(d.keys()):
                d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_COMPHET'].sum()]
                d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')) == w.replace('*','P')+'_NOCOMPHET'].sum()]

            otheralleleslength = len(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:])
            otherallelelist = list(rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['alleles'].values)
            if otheralleleslength > 1:
                #d['Other alleles ' + '(N=' + str(otheralleleslength) + ')'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                d['Other alleles'] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                d_comphet['Other alleles_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet['Other alleles_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
            else:
                theother = rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['alleles'].values[0]
                d[theother] = [rf.sort_values(by=['count'], ascending=False).reset_index().iloc[2:]['count'].sum()]
                d_comphet[theother +'_COMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_COMPHET' for i in otherallelelist]))].sum()]
                d_comphet[theother +'_NOCOMPHET'] = [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','P')).str.contains('|'.join([i.replace('*','P')+'_NOCOMPHET' for i in otherallelelist]))].sum()]
            genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
            genecountdf_comhet = pd.DataFrame(d_comphet)
            #genecountdf['from'] = country
        elif len(posalleles) == 2:
            d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[[0,1]]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[[0,1]]['count'])]))
            genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
            for w in list(d.keys()):
                d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','.')) == w.replace('*','.')+'_COMPHET'].sum()]
                d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','.')) == w.replace('*','.')+'_NOCOMPHET'].sum()]
            genecountdf_comhet = pd.DataFrame(d_comphet)
            #genecountdf['from'] = country
        elif len(posalleles) == 1:
            d = dict(zip(list(rf.sort_values(by=['count'], ascending=False).iloc[[0]]['alleles']),[[j] for j in list(rf.sort_values(by=['count'], ascending=False).iloc[[0]]['count'])]))
            for w in list(d.keys()):
                d_comphet[w+'_COMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','.')) == w.replace('*','.')+'_COMPHET'].sum()]
                d_comphet[w+'_NOCOMPHET']= [rf_comphet['count'].loc[rf_comphet['alleles'].apply(lambda x: x.replace('*','.')) == w.replace('*','.')+'_NOCOMPHET'].sum()]
            genecountdf = pd.DataFrame(d) #({'alleles':allele_names,GENE:allele_counts})
            genecountdf_comhet = pd.DataFrame(d_comphet)
            #genecountdf['from'] = country
    #rf = pd.DataFrame({'Gene':[GENE]*len(allele_names),'alleles':allele_names,'count':allele_counts})
    rf_aux = pd.DataFrame({'Gene':[GENE]*len(allele_names),'alleles':allele_names,'count':allele_counts_aux, 'country': [country]*len(allele_counts_aux), 'nr_found': [nr_found]*len(allele_counts_aux)})
    return genecountdf, genecountdf_comhet, otherallelelist, rf_aux


##
COLORS = ['tab:blue','tab:red','tab:green', 'tab:olive', 'tab:cyan', 'tab:orange', 'tab:brown','m','tab:pink']
othercountries = list(set(df['from']))
othercountries.remove('Spain')
othercountries.remove('Colombia')
othercountries.remove('Brazil')
othercountries = '|'.join(othercountries)
allcountries = '|'.join(list(set(df['from'])))
countries = [allcountries, 'Spain', 'Colombia', 'Brazil']#, othercountries]

RF = pd.DataFrame({'Gene':['hey'],'alleles':['what'], 'count':[0], 'country': [0], 'nr_found': [0]}) # initiallising this dataframe
fig, ax1 = plt.subplots(nrows, ncols, figsize=(7,10), sharey=False, sharex=False)

for GENE in GENES:
    genecountdf = pd.DataFrame()
    genecountdf_comhet = pd.DataFrame()
    coord_row = d_coords[GENE][1]
    coord_col = d_coords[GENE][0]
    country_col = list()
    otherallelistpergene = list()
    print(GENE)
    for country in countries:
        dff = df.loc[df['from'].str.contains(country)].groupby(GENE)['sample'].count().reset_index()
        dff = dff.loc[~(dff[GENE] == '')].copy()
        if len(dff) == 0:
            genecountdf = genecountdf.append(pd.Series(name=0))
            genecountdf_comhet = genecountdf_comhet.append(pd.Series(name=0))
            countryname = country.replace(allcountries, 'All').replace(othercountries, 'Other')
            country_col.append(countryname)
        else:
            if GENE == 'G6PD':
                genecountdf_aux, genecountdf_comhet_aux, otherallelelist_aux, rf = prepareDF4plot_chrX(df, country, GENE, RF)
            else:
                genecountdf_aux, genecountdf_comhet_aux, otherallelelist_aux, rf = prepareDF4plot(df, country, GENE, RF)
            otherallelistpergene.append(otherallelelist_aux)
            RF = pd.concat([RF,rf])
            countryname = country.replace(allcountries, 'All').replace(othercountries, 'Other')
            country_col.append(countryname)
            genecountdf = pd.concat([genecountdf, genecountdf_aux])
            genecountdf_comhet = pd.concat([genecountdf_comhet, genecountdf_comhet_aux])
        #print(country + '  ' + ' '.join(otherallelelist_aux))
    country_colcount = list()
    if coord_row == nrows-1:
        for c, cc in zip(country_col,countries):
            c = c + '\n'
            #country_colcount.append(c + '(n=' + str(len(df.loc[(~(df[GENE]).isnull()) & (df['from'].str.contains(cc))])) + ')')
            country_colcount.append(c + '[' + str(len(df.loc[(~(df[GENE]).isnull()) & (df['from'].str.contains(cc))]))+ ']')
    else:
        for c, cc in zip(country_col,countries):
            c = ' '
            #country_colcount.append(c + '(n=' + str(len(df.loc[(~(df[GENE]).isnull()) & (df['from'].str.contains(cc))])) + ')')
            country_colcount.append(c + '[' + str(len(df.loc[(~(df[GENE]).isnull()) & (df['from'].str.contains(cc))]))+ ']')
    genecountdf['from'] = country_colcount
    genecountdf_comhet['from'] = country_colcount
    nr_individuals = len(df.loc[~(df[GENE]).isnull()])
    f = list()
    othercheck = False
    for i in list(genecountdf.columns):
        if (i != 'Other alleles') and (i != 'from'):
            f.append(i)
        elif i == 'Other alleles':
            othercheck = True
    colors_aux = COLORS[0:len(f)]
    if othercheck:
        f.append('Other alleles')
        colors_aux.append('tab:purple')
    f_comhet = [[i+'_COMPHET',i+'_NOCOMPHET'] for i in f]
    f_comhet = [item for sublist in f_comhet for item in sublist]
    f.append('from')
    f_comhet.append('from')
    colors = colors_aux
    genecountdf = genecountdf[f].copy()
    genecountdf_comhet = genecountdf_comhet[f_comhet].copy()
    ylim = [0,100]
    otherallelistpergene_squeezed= list(set([item for sublist in otherallelistpergene for item in sublist]))
    genecountdf = genecountdf.rename(columns={'Other alleles':'Other [' + str(len(otherallelistpergene_squeezed)) + ']'})
    genecountdf.plot(x='from',ax = ax1[coord_row,coord_col],align = 'edge', kind = 'bar',width = -0.35, stacked = True, edgecolor = 'black', title = GENE , mark_right = True, color=colors)#,position=1+ ' (N=' + str(nr_individuals) + ')'
    label_list = list()
    for t in ax1[coord_row,coord_col].get_legend_handles_labels():
        label_list.append(t)
    #ax2 = ax1[coord_row,coord_col].twinx() # Create another axes that shares the same x-axis as ax.
    #genecountdf_comhet.plot(x='from',ax = ax2,position=0 ,kind = 'bar',width = 0.35, stacked = True, edgecolor = 'black' , mark_right = True, color=['lightgray','white']*len(colors)) #, title = GENE
    genecountdf = genecountdf.rename(columns={'Other alleles_COMPHET':'Other [' + str(len(otherallelistpergene_squeezed)) + ']','Other alleles_NOCOMPHET':'Other [' + str(len(otherallelistpergene_squeezed)) + ']'})
    genecountdf_comhet.plot(x='from',ax =  ax1[coord_row,coord_col] , align= 'edge',kind = 'bar',width = 0.15, stacked = True, edgecolor = 'black' , mark_right = True, color=['lightgray','white']*len(colors)) #, title = GENE
    if GENE == 'UGT1A1': # Adding a second legend
        label_list_comphet = [ax1[coord_row,coord_col].get_legend_handles_labels()[0][-2:], ['Comp.Het.\nHomozygous', 'No Comp.Het.\nNo Homoz.']]#  ax2.get_legend_handles_labels()[1][0:2]]
        legend2 = plt.legend(handles=label_list_comphet[0], labels=label_list_comphet[1],loc='center right', fontsize = 9, labelspacing=0.15, handletextpad=0.2,handlelength=1,bbox_to_anchor=(0.0,-0.1))
        ax1[coord_row,coord_col].add_artist(legend2)
    label_list = [label_list[0],label_list[1]]
    #label_list = [label_list[0] + label_list_comphet[0],label_list[1] + label_list_comphet[1]]
    #ax1[coord_row,coord_col].legend().set_visible(False)
    ax1[coord_row,coord_col].legend(handles=label_list[0], labels=label_list[1],loc='center left', fontsize = 9, labelspacing=0.1, handletextpad=0.2, handlelength=1,bbox_to_anchor=(1.0,0.5))
    xlim = [-0.5,3.3]
    ax1[coord_row, coord_col].set_xlim(xlim)
    ax1[coord_row, coord_col].set_ylim(ylim)
    ax1[coord_row, coord_col].set_title(label = GENE , fontsize= 10, fontweight="bold") #+ ' (N=' + str(nr_individuals) + ')', fontsize= 10) # title of plot
    plt.setp(ax1[coord_row, coord_col].get_xticklabels(), visible=True, rotation=0, ha='center')
    ax1[coord_row, coord_col].tick_params(axis = 'x',which='major', labelsize=10) #,labelrotation=0
    ax1[coord_row, coord_col].set_xlabel('')
    if coord_col == 0:
        ax1[coord_row,coord_col].set_ylabel('%')
    else:
        ax1[coord_row,coord_col].yaxis.set_visible(False)
plt.subplots_adjust(left=0.080, bottom=0.050, right=0.850, top=0.97, wspace=0.665, hspace=0.365)
plt.savefig('/path/to/Figures/Figure_2B_splittedbyLatamCountry.png',format = 'png', dpi = 500)
plt.show()
