#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:27:02 2020

@author: ckadelka
"""

##load libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import hierarchial_clustering as hc
from matplotlib import cm
import os
import itertools
import pareto

from metapub import PubMedFetcher #https://pypi.org/project/metapub/

filename = 'keywords_for_pubmed_search06'

EXACT_KEYWORD = True
GENERATE_JACCARD_PLOTS = False
def jaccard_similarity(x,y):
    intersection_cardinality = len(set.intersection(*[set(x), set(y)]))
    union_cardinality = len(set.union(*[set(x), set(y)]))
    return intersection_cardinality/float(union_cardinality)

#PMID of models from cell collective 
pmids_cellcollective = []
for fname in os.listdir('update_rules_cell_collective/'):
    if fname.endswith('.txt'):
        pmids_cellcollective.append(fname.split('_')[-1].split('.')[0])

#load the list of keywords from csv file
keywords = pd.read_csv(filename+'.csv',header=None)
keywords = np.array(keywords.iloc[:,0])
n_keywords = len(keywords)

#query Pubmed to find all articles (PMID) that contain the keyword
try:
    fetch = PubMedFetcher(api_key='7188286b0625041fa9c4be789db3665d3f08')
except TypeError:
    fetch = PubMedFetcher()


pmids_per_keyword = []
category1 = ['Boolean','logical','dynamical','dynamic','logic']
category2 = ['regulatory','signalling','signaling','transduction','network','networks']
category3 = ['network','networks','circuit','model','models','modelled','modeled','modelling','model','pathway','pathways']
   
category1 = ['Boolean','logical','dynamical','dynamic','logic']
category2 = ['regulatory','signalling','signaling','network','networks']
category3 = ['network','networks','circuit','model','models','modelled','modeled','modelling','model','pathway','pathways']
    
 
keywords_categories = []
pmids_per_keyword_combination = []
for key1 in category1:
    for key2 in category2:
        for key3 in category3:
            if key2 == key3: #network/networks
                keyword_combination = key1+' '+key2
            elif key2 in key3 or key3 in key2:
                continue
            else:
                keyword_combination = key1+' '+key2+' '+key3
            keywords_categories.append(keyword_combination)
            articles = fetch.pmids_for_query(keyword_combination,retstart=0, retmax=1000000)
            #articles = set(pmids_per_keyword[keywords.index(key1)]) & set(pmids_per_keyword[keywords.index(key2)]) & set(pmids_per_keyword[keywords.index(key3)])
            pmids_per_keyword_combination.append(articles)
        print(key2)
    print(key1)
        
    
pairs =['"%s"' % el for el in '''boolean circuit
synchronous boolean
boolean model
logical modelling
boolean modelling
boolean modeling
logical modeling
apoptosis network
perturbation simulations
input combinations
logical model
boolean network
network describing
dynamic boolean model
patterns observed experimentally
boolean network model
underlying regulatory network
regulatory network model
budding yeast cell
regulatory network controlling'''.split('\n')]

for word in pairs:
    keywords_categories.append(word)
    articles = fetch.pmids_for_query(word,retstart=0, retmax=1000000)
    #articles = set(pmids_per_keyword[keywords.index(key1)]) & set(pmids_per_keyword[keywords.index(key2)]) & set(pmids_per_keyword[keywords.index(key3)])
    pmids_per_keyword_combination.append(articles)

n_pmids_per_keyword_combination = list(map(len,pmids_per_keyword_combination))

if GENERATE_JACCARD_PLOTS:
    #look at the similarity of the results based on the different keywords
    sim_jaccards = np.ones((n_keywords,n_keywords))
    for i in range(n_keywords):
        for j in range(i+1,n_keywords):
            sim_jaccards[i,j] = jaccard_similarity(pmids_per_keyword[i],pmids_per_keyword[j])
            sim_jaccards[j,i] = sim_jaccards[i,j]
    
    labels = ['%s\n(%i)' % (a,b) for a,b in zip(keywords,n_pmids_per_keyword)]
    hc.heatmap(sim_jaccards,labels,labels,'average','average','euclidean','euclidean',cm.Greens,'sim_jaccard_%s.pdf' % (filename+('_exact' if EXACT_KEYWORD else '')),'Jaccard Similarity',CENTERED_AT_0=False)

#how many models from cell collective are found by each keyword
pmids_overlap_per_keyword = np.array([set(el) & set(pmids_cellcollective) for el in pmids_per_keyword_combination])
n_pmids_overlap_per_keyword = list(map(len,pmids_overlap_per_keyword))

#how many models from cell collective are found total
all_retrieved_pmids = set.union(*[set([])] + [set(el) for el in pmids_per_keyword_combination])
pmids_overlap_all_keywords = set(all_retrieved_pmids) & set(pmids_cellcollective)
n_pmids_overlap_all_keywords = len(pmids_overlap_all_keywords)

#display those articles from cell collective that we missed
missed_articles = []
for PMID in set(pmids_cellcollective)-pmids_overlap_all_keywords:
    article = fetch.article_by_pmid(PMID)
    missed_articles.append(article)
    if 1>0:#'Boolean' in article.title+'\n'+article.abstract:
        print((article.title))
        print((article.abstract))
        print()
        
#Jaccard index per keyword combo
res = []
n_CC = len(set(pmids_cellcollective))
jaccards = []
ratios = []
for i in range(len(n_pmids_overlap_per_keyword)):
    count = n_pmids_overlap_per_keyword[i]
    n_articles = n_pmids_per_keyword_combination[i]
    if n_articles==0:
        continue
    word = keywords_categories[i]
    ratios.append(count/n_articles)
    jaccards.append(count/(n_articles+n_CC-count))
    res.append([word,count,n_articles,count/n_articles,count/(n_articles+n_CC-count)])
A = pd.DataFrame(res,columns=['word','in CC','in Pubmed','ratio','JI'])
A.sort_values('JI',ascending=False).to_excel('occurence_of_word_combos.xlsx')


dict_score = dict()
dict_keyword = dict()
for i in range(len(n_pmids_overlap_per_keyword)):
    count = n_pmids_overlap_per_keyword[i]
    n_articles = n_pmids_per_keyword_combination[i]
    word = keywords_categories[i] 
    for article in pmids_per_keyword_combination[i]:
        
        try:
            score = dict_score[article]
            if score<=count/n_articles:
                dict_keyword[article] = word
                dict_score[article] = count/n_articles
        except KeyError:
            dict_score.update({article:count/n_articles})
            dict_keyword.update({article:word})

models_found_2020_03_02 = [26244885,16873462,18250321,22267503,17722974,20221256,20862356,21563979,19118495,25538703,22871178,18852469,22102804,26102287,22253585,18433497,23171249,26446703,18463633,21968890,26090929,23743337,24970389,25431332,26163548,26207376,26340681,26385365,26528548,26573569,26616283,28639170,29513758,17010384,19144179,19185585,22962472,23056457,23134720,23233838,23520449,23868318,24250280,25908096,26088082,27148350,27542373,27594840,28361666,29206223,16464248,19662154,16542429,19025648,19422837,23049686,23081726,23637902,25790483,26408858,26751566,28584084,30518777,26660865,22448278,27765040,28381275,16968132,19524598,23169817,30024932,30953496,1753781,11082279,30038409,30530226,24564942,26603105,29230182,18194572,23935937,24069138,24079299,25780058,25980672,27464342,27613445,28426669,20169167,30104572,18301750,20659480,26346668,29237040,29622038,21639591,23764028,23469179,22979979,20011108,20824124,24550717,24729835,15486106,19763337,20655356,21342546,26259567,29107989,12957124,18689831,18692073,18956339,23463387,23922675,25398016,27229461,27587700,27599298,28155685,28303788,28968438,30053801,30116200,30546316,25163068,12782112,29671404,18978941,23017186,23425857,23813012,25768678,28663244]
models_found_2020_03_02 = list(map(str,models_found_2020_03_02))
PMIDS_investigated_2020_03_02 = [26244885,1753781,11082279,12782112,12957124,15486106,16464248,16542429,16873462,16968132,17010384,17722974,18194572,18250321,18301750,18433497,18463633,18689831,18692073,18852469,18956339,18978941,19025648,19118495,19144179,19185585,19422837,19524598,19662154,19763337,20011108,20169167,20221256,20655356,20659480,20824124,20862356,21342546,21563979,21639591,21968890,22102804,22253585,22267503,22448278,22871178,22962472,22979979,23017186,23049686,23056457,23081726,23134720,23169817,23171249,23233838,23425857,23463387,23469179,23520449,23637902,23743337,23764028,23813012,23868318,23922675,23935937,24069138,24079299,24250280,24550717,24564942,24729835,24970389,25163068,25398016,25431332,25538703,25768678,25780058,25790483,25908096,25980672,26088082,26090929,26102287,26163548,26207376,26259567,26340681,26346668,26385365,26408858,26446703,26528548,26573569,26603105,26616283,26660865,26751566,27148350,27229461,27464342,27542373,27587700,27594840,27599298,27613445,27765040,28155685,28303788,28361666,28381275,28426669,28584084,28639170,28663244,28968438,29107989,29206223,29230182,29237040,29513758,29622038,29671404,30024932,30038409,30053801,30104572,30116200,30518777,30530226,30546316,30953496,10223,15879,435523,469416,528944,623579,636866,711992,851664,903872,986544,1102721,1238113,1268282,1339957,1344673,1387400,1401550,1402706,1443287,1497753,1560184,1573904,1603808,1623291,1697772,1736362,1742982,1747785,1775661,1793756,1827831,2175124,2663334,2758786,2974865,3267769,3950540,5297650,5883205,6227688,6510054,6605796,6997678,7358922,7550157,7673867,7673868,7790010,7844509,7885260,7889762,7968557,8027704,8294127,8332424,8440975,8447506,8525870,8609299,8653099,8693955,8805827,8817903,8920841,8961824,9090153,9398508,9413865,9488886,9506386,9801347,10011547,10014526,10020383,10025541,10032277,10039533,10046267,10053243,10056406,10059380,10063130,10146119,10152882,10169569,10194925,10194926,10205599,10235650,10335601,10373357,10380182,10380189,10388754,10394139,10404976,10431185,10438671,10439869,10441458,10456665,10460674,10485897,10487867,10501921,10521334,10528071,10562736,10615501,10636026,10636093,10643753,10655182,10668659,10675717,10697599,10704301,10707821,10715234,10732490,10733855,10743952,10751456,10790766,10810010,10824422,10824431,10828326,10837909,10866829,10902178,10902179,10940402,10950837,10955645,10977069,10979191,10988040,10991019,10995349,11008375,11015189,11015957,11019016,11019303,11019360,11021655,11031615,11042046,11042050,11072331,11073207,11079959,11088116,11088285,11088334,11088448,11088804,11089619,11099257,11099258,11106956,11108466,11112254,11114226,11138149,11138166,11171388,11181111,11182199,11193103,11233025,11233511,11261316,11262079,11262963,11267737,11267738,11270747,11290142,11291730,11291941,11303522,11304372,11304386,11304390,11308646,11308814,11311467,11315315,11318171,11328680,11336120,11340853,11341566,11343406,11349518,11371430,11414939,11415183,11415228,11461228,11497549,11497645,11516449,11518665,11519450,11535110,11539411,11547977,11553811,11569498,11583529,11586365,11595328,11668442,11679699,11681752,11681759,11690098,11690138,11697671,11700600,11708322,11714557,11716864,11734025,11735289,11736320,11754340,11766939,11786029,11800715,11800742,11800758,11847074,11911796,11934021,11942532,11970520,12005947,12049761,12135461,12241409,12715183,12804224,12906408,12920140,12923549,12934604,12935020,12940736,12946201,12957633,12963822,14604408,14615204,14657375,14724639,14756587,15068191,15189048,15244911,15246788,15478442,15495907,15520543,15572773,15736517,15742883,15783412,15820736,15903494,15906502,16132344,16153042,16485968,16846971,16870321,16986551,16986617,17069276,17109765,17155099,17221085,17316697,17543998,17552086,17646673,17653397,17945690,17946091,18255652,18267419,18269749,18321113,18359043,18384716,18488668,18517473,18550220,18614585,18626127,18791235,18792497,19014577,19014957,19039389,19113843,19254727,19347209,19351582,19426782,19477975,19499341,19517238,19584836,19588107,19649409,19658525,19684935,19732446,19750394,19785753,19792032,19884180,19905381,20004671,20014476,20024808,20126540,20222348,20232996,20329858,20364976,20478915,20531402,20601441,20628620,20646315,20665073,20704717,20707925,20808719,20863397,20885776,20953502,21121062,21198122,21212881,21214342,21230358,21311809,21325303,21417939,21464507,21489274,21554763,21685052,21702798,21808895,21828085,21829342,21918620,21977695,22016404,22128010,22132067,22150915,22151781,22166257,22174668,22212351,22239642,22300322,22303372,22356156,22443451,22534843,22556185,22558096,22584015,22641521,22645556,22729247,22849591,22879903,22927416,22929591,22932419,22938561,22941654,22953690,23011283,23020215,23027966,23033959,23034351,23096153,23136077,23142247,23277042,23320651,23335016,23361989,23565141,23670445,23679474,23761998,23822502,23822512,23835289,23975240,24023735,24028533,24047877,24130232,24176122,24194720,24212100,24213870,24339029,24351986,24388494,24391592,24444824,24451547,24484917,24631047,24651476,24884358,24886608,24932713,24983623,25078949,25080304,25093019,25116096,25122140,25184579,25258177,25269159,25315877,25320069,25324091,25345703,25433558,25462327,25482233,25493491,25493675,25573782,25684199,25705700,25723815,25768836,25796455,25819791,25847279,25909953,25954815,25967891,25972988,25994715,26061094,26078349,26106462,26111433,26132266,26231087,26238783,26280446,26290194,26305513,26390498,26428008,26439385,26452000,26465843,26496494,26671814,26701883,26818802,26851011,26857675,26981147,27105325,27138333,27306057,27413618,27417985,27504165,27587699,27600248,27635462,27716031,27760125,27774993,27932399,27932972,27965080,27978857,27993914,28073755,28113565,28150999,28155725,28178334,28182683,28265989,28268604,28269650,28362589,28402858,28425119,28513208,28548604,28584110,28614451,28640804,28649434,28663080,28724002,28726971,28736575,28792568,28814968,28881959,28885158,28904202,28953896,29038942,29133455,29158122,29178844,29212542,29267677,29267816,29267823,29322926,29393068,29428070,29476134,29562019,29566771,29670664,29671395,29684020,29705481,29740342,29758614,29784817,29792308,29872402,29875674,29879983,29924153,29937735,29949973,29971008,29988359,29989991,30001194,30034343,30036529,30169736,30306904,30319440,30384638,30384642,30415697,30421414,30606121,725138,1824522,6241491,22081601,26535546]
PMIDS_investigated_2020_03_02 = list(map(str,PMIDS_investigated_2020_03_02))

all_articles = list(dict_score.keys())
all_articles_position = dict(zip(all_articles,list(range(len(all_articles)))))
all_articles_in_CC = [0 for article in all_articles]
for PMID in pmids_cellcollective:
    try:
        all_articles_in_CC[all_articles_position[PMID]] = 1
    except KeyError: #article from CC not found
        pass
all_articles_scores = [dict_score[article] for article in all_articles]
all_articles_words = [dict_keyword[article] for article in all_articles]
all_articles_contains_model = [0 for article in all_articles]
for PMID in models_found_2020_03_02:
    try:
        all_articles_contains_model[all_articles_position[PMID]] = 1
    except KeyError: #article from CC not found
        pass
all_articles_already_investigated = [0 for article in all_articles]
for PMID in PMIDS_investigated_2020_03_02:
    try:
        all_articles_already_investigated[all_articles_position[PMID]] = 1
    except KeyError: #article from CC not found
        pass    
 
#export option 1:
score_per_article = pd.DataFrame(np.c_[all_articles,all_articles_scores,all_articles_words,all_articles_in_CC,all_articles_contains_model,all_articles_already_investigated],columns=['PMID','max ratio','word','CC','model','investigated'])
score_per_article.sort_values('max ratio',ascending=False).to_excel('occurence_of_word_combos_scoring.xlsx')


#export option 2:
A = pd.read_excel('Pubmed_GRN_Database_Progress_sheet.xlsx')
count = len(A)
A['score'] = pd.Series(['' for el in range(count)], index=A.index)
A['best keyword'] = pd.Series(['' for el in range(count)], index=A.index)
n_columns = len(A.columns)
for i in range(len(A)):
    try:
        index = all_articles_position[str(A['PMID'][i])]
        A.loc[i,'score'] = all_articles_scores[index]
        A.loc[i,'best keyword'] = all_articles_words[index]        
    except KeyError: #PMID in excel file is not among best IDs
        A.loc[i,'score'] = 0
        pass

#A = A.sort_values('PMID',ascending=True)
#A = A.sort_values('No, does not',ascending=False)
#A = A.sort_values('Unsure',ascending=False)
A = A.sort_values(['Yes, contains model','Unsure','No, does not','PMID'],ascending=[False,False,False,True])

dict_id_in_excel = dict(zip(np.array(A['PMID']),range(len(A))))
for i in range(len(all_articles)):
    try:
        dict_id_in_excel[int(all_articles[i])]
    except KeyError:
        index = count
        A.loc[index] = [np.int(all_articles[i])] + [np.nan for el in range(n_columns-1)]
        A.loc[index,'score'] = all_articles_scores[i]
        A.loc[i,'best keyword'] = all_articles_words[i]
        count+=1
    if i%500==0:
        print(i)
    if i>5000:
        break
        
A.to_excel('Pubmed_GRN_Database_Progress_sheet_mod.xlsx')


#export option 2b (as np-array):
A = pd.read_excel('Pubmed_GRN_Database_Progress_sheet.xlsx')
A_np = np.array(A)
count = A_np.shape[0]
A_np = np.c_[A_np,np.nan*np.ones((count,2))]
n_columns = A_np.shape[1]
pmids = A_np[:,list(A.columns).index('PMID')].reshape(count)
modelyes = A_np[:,list(A.columns).index('Yes, contains model')].reshape(count)
modelmaybe = A_np[:,list(A.columns).index('Unsure')].reshape(count)
modelno = A_np[:,list(A.columns).index('No, does not')].reshape(count)
modelyes[np.bitwise_and(modelyes!='x',modelyes!='X')]=''
modelno[np.bitwise_and(modelno!='x',modelno!='X')]=''
modelmaybe[np.bitwise_and(modelmaybe!='x',modelmaybe!='X')]=''

scores = []
best_keywords = []
for i in range(len(A)):
    try:
        index = all_articles_position[str(pmids[i])]
        scores.append(all_articles_scores[index])
        best_keywords.append(all_articles_words[index])        
    except KeyError: #PMID in excel file is not among best IDs
        scores.append(0)
        best_keywords.append('')  


A_np[:,-2] = scores
A_np[:,-1] = best_keywords
indices = np.array(sorted(list(range(count)),key=lambda x: (modelyes[x],modelmaybe[x],modelno[x],1e13-np.array(pmids[x],dtype=int)),reverse=True))
#scores = np.array(scores)[indices]
#best_keywords = np.array(best_keywords)[indices]
pmids = pmids[indices]
A_np = A_np[indices,:]

dict_id_in_excel = dict(zip(pmids,range(count)))
B_np = []
for i in range(len(all_articles)):
    try:
        dict_id_in_excel[int(all_articles[i])]
    except KeyError:
        B_np.append([(all_articles[i])]+['' for el in range(A_np.shape[1]-3)]+[all_articles_scores[i],all_articles_words[i]])
    if i%500==0:
        print(i)
    if i>5000:
        break
B_np = np.array(B_np)
C_np = np.r_[A_np,np.array(B_np)]

pd.DataFrame(C_np,columns=list(A.columns)+['score','best keyword']).to_excel('Pubmed_GRN_Database_Progress_sheet_mod2.xlsx')





#best pairs/triples/etc of keywords, i.e. lowest number of false positives / PMIDs
keyword_combinations = []
n_pmids_per_keyword_combination = []
n_pmids_overlap_per_keyword_combination = []
for no_keywords in range(1,11):
    for keyword_ids in itertools.combinations(list(range(n_keywords)), no_keywords):
        keyword_combinations.append(keyword_ids)
        dummy = set.union(*[set([])] + [set(pmids_per_keyword[i]) for i in keyword_ids])
        n_pmids_per_keyword_combination.append( len(dummy) )
        n_pmids_overlap_per_keyword_combination.append( len(set(dummy)&set(pmids_cellcollective)) )
        
f,ax = plt.subplots()
ax.semilogx(n_pmids_per_keyword_combination,n_pmids_overlap_per_keyword_combination,'o')
ax.set_xlabel('No PMIDs sharing at least one keyword')
ax.set_ylabel('No PMIDs recovered from Cell Collective')
pareto_frontier = pareto.get_pareto_frontier(np.c_[-np.array(n_pmids_per_keyword_combination),n_pmids_overlap_per_keyword_combination])
pareto_frontier[:,0] *= -1
ax.semilogx(pareto_frontier[:,0],pareto_frontier[:,1],'ro')
plt.savefig('No_PMIDs_recovered_from_cellcollective_%s.pdf' % (filename+('_exact' if EXACT_KEYWORD else '')))

dummy = list((a,b,c) for a,b,c in zip(keyword_combinations,n_pmids_per_keyword_combination,-np.array(n_pmids_overlap_per_keyword_combination)))
dummy.sort( key = lambda dummy: (dummy[2],dummy[1]) )

pareto_optimal = [dummy[0][0]]
pareto_optimal_x = [dummy[0][1]]
pareto_optimal_y = [-dummy[0][2]]
pareto_optimal_PMIDs = [set.union(*[set([])] + [set(pmids_per_keyword[ii]) for ii in dummy[0][0]])]
currrent_val = -1
best_x_thus_far = 1e10
for i in range(1,len(dummy)):
    if dummy[i][2]>dummy[i-1][2] and best_x_thus_far>dummy[i][1]:
        best_x_thus_far = dummy[i][1]
        pareto_optimal.append(dummy[i][0])
        pareto_optimal_x.append(dummy[i][1])
        pareto_optimal_y.append(-dummy[i][2])
        pareto_optimal_PMIDs.append(set.union(*[set([])] + [set(pmids_per_keyword[ii]) for ii in dummy[i][0]]))
    else:
        continue

OUT = pd.DataFrame(np.c_[pareto_optimal_y,pareto_optimal_x,[' OR '.join([keywords[kk] for kk in k]) for k in pareto_optimal]],columns = ['No PMIDs recovered from Cell Collective','No PMIDs sharing at least one keyword','keyword combination'])
OUT.to_excel('pareto_efficient_search_terms_%s.xlsx' % (filename+('_exact' if EXACT_KEYWORD else '')))

#look at all titles and abstracts in cell collective and collect common words
f = open('most_common_3000_english_words.txt','r')
common_words = f.read().splitlines()
f.close()
common_words.extend(['is','has','was','been','an','are','were'])
common_words.remove('model')
common_words.remove('network')
common_words.remove('cell')
common_words.remove('gene')
common_words_dict = dict(list(zip(common_words,list(range(len(common_words))))))

uncommon_words = []
uncommon_words_dict = {}
count_uncommon_words = []
n_uncommon_words = 0
BN = 0
for PMID in set(pmids_cellcollective):
    dummy_count_uncommon_words = [0 for _ in range(n_uncommon_words)]
    try:
        f = open('title_abstract_cell_collective/title_abstract_%s.txt' % PMID, 'r')
        text = f.read()#.decode('utf8')
        f.close()
    except IOError: #if file hasn't been created
        article = fetch.article_by_pmid(PMID)
        f = open('title_abstract_cell_collective/title_abstract_%s.txt' % PMID, 'w')
        text = article.title + '\n' + article.abstract
        f.write(text)#.encode('utf8'))
        f.close()
        
    BN += int('Boolean' in text and 'network' in text)
    words = text.replace('.','').replace(';','').replace(',','').replace(':','').replace('-','').replace('"','').split(' ')
    for word in words:
        word = word.lower()
        try: #check whether the word is common
            common_words_dict[word]
        except KeyError:
            try: #check whether the word is uncommon but already been found before
                index = uncommon_words_dict[word]
                dummy_count_uncommon_words[index] = 1
            except: #new uncommon word
                uncommon_words_dict.update( {word : n_uncommon_words} )
                uncommon_words.append(word)
                count_uncommon_words.append(0)
                dummy_count_uncommon_words.append(1)
                n_uncommon_words+=1
    count_uncommon_words = [a+b for a,b in zip(count_uncommon_words,dummy_count_uncommon_words)]

dummy = list((a,b) for a,b in zip(uncommon_words,-np.array(count_uncommon_words)))
dummy.sort( key = lambda dummy: (dummy[1]) )

res = []
min_number = 7
n_CC = len(set(pmids_cellcollective))
for (word,count) in dummy:
    if -count<=min_number:
        break
    else:
        articles = fetch.pmids_for_query('"%s"' % word,retstart=0, retmax=1000000)
        n_articles = len(articles)
        print(word,-count,n_articles,-count/n_articles,-count/(n_articles+n_CC+count))
        res.append([word,-count,n_articles,-count/n_articles,-count/(n_articles+n_CC+count)])
A = pd.DataFrame(res,columns=['word','in CC','in Pubmed','ratio','JI'])
A.sort_values('JI',ascending=False).to_excel('occurence_of_words_CC_vs_Pubmed_at_least_%i_times.xlsx' % min_number)
    




#pairs of words
uncommon_words = []
uncommon_words_dict = {}
count_uncommon_words = []
n_uncommon_words = 0
BN = 0
for PMID in set(pmids_cellcollective):
    dummy_count_uncommon_words = [0 for _ in range(n_uncommon_words)]
    try:
        f = open('title_abstract_cell_collective/title_abstract_%s.txt' % PMID, 'r')
        text = f.read()#.decode('utf8')
        f.close()
    except IOError: #if file hasn't been created
        article = fetch.article_by_pmid(PMID)
        f = open('title_abstract_cell_collective/title_abstract_%s.txt' % PMID, 'w')
        text = article.title + '\n' + article.abstract
        f.write(text)#.encode('utf8'))
        f.close()
    
    words = np.array(list(map(str.lower,text.replace('.',' endsentence ').replace(';','').replace('\n','').replace(',','').replace(':','').replace('-','').replace('"','').split(' '))))
    words = words[words!='']
    for ii in range(len(words)-1):
        if words[ii]=='endsentence' or words[ii+1] == 'endsentence':
            continue
        counter_common_words = 0
        try:
            common_words_dict[words[ii]]
            counter_common_words+=1
        except KeyError:
            pass
        try:
            common_words_dict[words[ii+1]]
            counter_common_words+=1
        except KeyError:
            pass
        if counter_common_words>0:
            continue
        wordpair = words[ii]+' '+words[ii+1]
        try: #check whether the word is uncommon but already been found before
            index = uncommon_words_dict[wordpair]
            dummy_count_uncommon_words[index] = 1
        except: #new uncommon word
            uncommon_words_dict.update( {wordpair : n_uncommon_words} )
            uncommon_words.append(wordpair)
            count_uncommon_words.append(0)
            dummy_count_uncommon_words.append(1)
            n_uncommon_words+=1
    count_uncommon_words = [a+b for a,b in zip(count_uncommon_words,dummy_count_uncommon_words)]

dummy = list((a,b) for a,b in zip(uncommon_words,-np.array(count_uncommon_words)))
dummy.sort( key = lambda dummy: (dummy[1]) )

res = []
min_number = 1
n_CC = len(set(pmids_cellcollective))
for (word,count) in dummy:
    if -count<=min_number:
        break
    else:
        articles = fetch.pmids_for_query('"%s"' % word,retstart=0, retmax=1000000)
        n_articles = len(articles)
        print(word,-count,n_articles,-count/n_articles,-count/(n_articles+n_CC+count))
        res.append([word,-count,n_articles,-count/n_articles,-count/(n_articles+n_CC+count)])
A = pd.DataFrame(res,columns=['word','in CC','in Pubmed','ratio','JI'])
A.sort_values('JI',ascending=False).to_excel('occurence_of_word_pairs_CC_vs_Pubmed_at_least_%i_times.xlsx' % (min_number+1))
    


#triplets of words
uncommon_words = []
uncommon_words_dict = {}
count_uncommon_words = []
n_uncommon_words = 0
BN = 0
for PMID in set(pmids_cellcollective):
    dummy_count_uncommon_words = [0 for _ in range(n_uncommon_words)]
    try:
        f = open('title_abstract_cell_collective/title_abstract_%s.txt' % PMID, 'r')
        text = f.read()#.decode('utf8')
        f.close()
    except IOError: #if file hasn't been created
        article = fetch.article_by_pmid(PMID)
        f = open('title_abstract_cell_collective/title_abstract_%s.txt' % PMID, 'w')
        text = article.title + '\n' + article.abstract
        f.write(text)#.encode('utf8'))
        f.close()
        
    words = np.array(list(map(str.lower,text.replace('.',' endsentence ').replace(';',' endsentence ').replace(',','').replace(' endsentence ','').replace('-','').replace('"','').replace('  ',' ').split(' '))))
    words = words[words!='']
    for ii in range(len(words)-2):
        if words[ii]=='endsentence' or words[ii+1] == 'endsentence' or words[ii+2]=='endsentence':
            continue
        counter_common_words = 0
        try:
            common_words_dict[words[ii]]
            counter_common_words+=1
        except KeyError:
            pass
        try:
            common_words_dict[words[ii+1]]
            counter_common_words+=1
        except KeyError:
            pass
        try:
            common_words_dict[words[ii+2]]
            counter_common_words+=1
        except KeyError:
            pass
        if counter_common_words>0:
            continue
        wordpair = words[ii]+' '+words[ii+1]+' '+words[ii+2]
        try: #check whether the word is uncommon but already been found before
            index = uncommon_words_dict[wordpair]
            dummy_count_uncommon_words[index] = 1
        except: #new uncommon word
            uncommon_words_dict.update( {wordpair : n_uncommon_words} )
            uncommon_words.append(wordpair)
            count_uncommon_words.append(0)
            dummy_count_uncommon_words.append(1)
            n_uncommon_words+=1
    count_uncommon_words = [a+b for a,b in zip(count_uncommon_words,dummy_count_uncommon_words)]

dummy = list((a,b) for a,b in zip(uncommon_words,-np.array(count_uncommon_words)))
dummy.sort( key = lambda dummy: (dummy[1]) )


res = []
min_number = 1
n_CC = len(set(pmids_cellcollective))
for (word,count) in dummy:
    if -count<=min_number:
        break
    else:
        articles = fetch.pmids_for_query('"%s"' % word,retstart=0, retmax=1000000)
        n_articles = len(articles)
        print(word,-count,n_articles,-count/n_articles,-count/(n_articles+n_CC+count))
        res.append([word,-count,n_articles,-count/n_articles,-count/(n_articles+n_CC+count)])
A = pd.DataFrame(res,columns=['word','in CC','in Pubmed','ratio','JI'])
A.sort_values('JI',ascending=False).to_excel('occurence_of_word_triplets_CC_vs_Pubmed_at_least_%i_times.xlsx' % (min_number+1))
    


#save abstracts to articles on the Pareto frontier
fetch = PubMedFetcher()
counter = 0
for pareto_optimal_PMID_list in pareto_optimal_PMIDs[::-1]:
    for PMID in pareto_optimal_PMID_list:
        try:
            try:
                f = open('title_abstract_pubmed/title_abstract_%s.txt' % PMID, 'r')
                f.close()
            except IOError: #if file hasn't been created
                article = fetch.article_by_pmid(PMID)
                f = open('title_abstract_pubmed/title_abstract_%s.txt' % PMID, 'w')
                text = article.title + '\n' + article.abstract
                f.write(text.encode('utf8'))
                f.close()
            counter +=1
            if counter%10 == 0:
                print(counter)
        except TypeError: #there are some PMIDs that link to papers without abstract (or possibly title), exclude them
            continue



f = open('progress_table.txt','w')
PMIDs_in_table_dict = {}
len_PMIDs_in_table_dict = 0
for ii,pareto_optimal_PMID_list in pareto_optimal_PMIDs[::-1]:
    for PMID in pareto_optimal_PMID_list:
        for kk,keyword in keywords:
            try:
                index = PMIDs_in_table_dict[PMID]
            except KeyError:
                PMIDs_in_table_dict.update({PMID:len_PMIDs_in_table_dict})















##New idea: only consider those abstracts where our keyword occurs in a sentence including an action verb like generated, created, etc
action_verb_parts = ['create','generate','present','propose','use','construct','translate','set','model','establish','provide','introduce','identify','build','compile','implement','develop','formalize','encode','reconstruct','integrate','postulate','derive','infer','convert','define']
action_verb_parts.extend( ['created','generated','presented','proposed','used','constructed','translated','set','modeled','established','provided','introduced','identified','built','compiled','implemented','developed','formalized','encoded','reconstructed','integrated','postulated','derived','inferred','converted','defined'])
action_verb_parts.extend( ['creating','generating','presenting','proposing','using','constructing','translating','setting','modeling','establishing','providing','introducing','identifying','building','compiling','implementing','developing','formalizing','encoding','reconstructing','integrating','postulating','deriving','inferring','converting','defining'])


fetch = PubMedFetcher()
counter = 0
pmids_per_keyword_extended = [el[:] for el in pmids_per_keyword] #deep copy, necessary for lists
n_keywords_found = np.zeros(n_keywords)
for kk,keyword in enumerate(keywords):
    pmids_per_keyword_extended.append([])
    for PMID in pmids_per_keyword[kk]:
        try:
            f = open('title_abstract_pubmed/title_abstract_%s.txt' % PMID, 'r')
            text = f.read().decode('utf8')
            f.close()
        except IOError: #if file hasn't been created
            article = fetch.article_by_pmid(PMID)
            f = open('title_abstract_pubmed/title_abstract_%s.txt' % PMID, 'w')
            text = article.title + '\n' + article.abstract
            f.write(text.encode('utf8'))
            f.close()
        except TypeError: #there are some PMIDs that link to papers without abstract (or possibly title), exclude them
            continue
        sentences = article.abstract.split('.')
        VERB_AND_KEYWORD = False
        FOUND_KEYWORD = False
        for sentence in sentences:
            if keyword in sentence:
                FOUND_KEYWORD = True
                if np.any([verb in sentence for verb in action_verb_parts]):
                    VERB_AND_KEYWORD = True
                    break
        if VERB_AND_KEYWORD:
            pmids_per_keyword_extended[-1].append(PMID)
        if FOUND_KEYWORD:
            n_keywords_found[kk] +=1





for keyword in keywords:
    article = fetch.article_by_pmid(PMID)

    articles = fetch.pmids_for_query(keyword,retstart=0, retmax=100000)
    for article in articles:
        sentences = articles.split('.')
    pmids_per_keyword.append(articles)
n_pmids_per_keyword = list(map(len,pmids_per_keyword))



