#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:01:32 2023

@author: ckadelka
"""

##Imports

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import canalizing_function_toolbox_v13 as can
import load_database13 as db

plt.rcParams.update({'font.size': 9})

# #to decide which models to keep for the meta-analysis we ran this and manually inspected the similar models regarding their overlap, variables, etc.
# folders=['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/']
# max_degree=2
# max_N=1000
# models_to_keep = []
# [Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded] = db.load_database(folders,max_degree=max_degree,max_N=max_N)
# similar_sets_jaccard = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=0.8,USE_JACCARD = True,models_to_keep=models_to_keep)[-1]
# similar_sets_overlap = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=0.9,USE_JACCARD = False,models_to_keep=models_to_keep)[-1]
# Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_excluded,similar_sets_overlap = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=0.9,USE_JACCARD = False,models_to_keep=models_to_keep)
# similar_sets_overlap = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=0.8,USE_JACCARD = False,models_to_keep=models_to_keep)[-1]

# variabless_simple = [[el.lower().replace('_','').replace('.','').replace('kappa','k').replace('-','') for el in variables] for variables in variabless]

# for cluster in similar_sets_overlap:
#     n = len(cluster)
#     indices = [models_loaded.index(el) for el in cluster]
#     print(indices,cluster)
#     for i in range(n):
#         for j in range(i+1,n):
#             sim_overlap = db.overlap_similarity(variabless_simple[indices[i]],variabless_simple[indices[j]])
#             sim_jaccard = db.jaccard_similarity(variabless_simple[indices[i]],variabless_simple[indices[j]])
#             vars_i_not_j = list(set(variabless_simple[indices[i]])-set(variabless_simple[indices[j]]))
#             vars_j_not_i = list(set(variabless_simple[indices[j]])-set(variabless_simple[indices[i]]))
#             vars_i_not_j.sort()
#             vars_j_not_i.sort()
#             print(i,j,len(variabless_simple[indices[i]]),len(variabless_simple[indices[j]]),len(constantss[indices[i]]),len(constantss[indices[j]]),sim_overlap,sim_jaccard,vars_i_not_j,vars_j_not_i)
#     print()

# pmids_models_loaded = [int(''.join([entry for entry in el.split('.')[0].split('_') if (entry[0] in '123456789' and len(entry)>3)])) for el in models_loaded]

# variabless_simple = [[el.lower().replace('_','').replace('.','').replace('kappa','k').replace('-','') for el in variables] for variables in variabless]


# a = pd.value_counts(pmids_models_loaded)
# for PMID in list(a.index)[:sum(a>1)]:
#     cluster=np.array(models_loaded)[can.find_all_indices(pmids_models_loaded,PMID)]
#     n = len(cluster)
#     indices = [models_loaded.index(el) for el in cluster]
#     print(indices,cluster)
#     for i in range(n):
#         for j in range(i+1,n):
#             sim_overlap = db.overlap_similarity(variabless_simple[indices[i]],variabless_simple[indices[j]])
#             sim_jaccard = db.jaccard_similarity(variabless_simple[indices[i]],variabless_simple[indices[j]])
#             vars_i_not_j = list(set(variabless_simple[indices[i]])-set(variabless_simple[indices[j]]))
#             vars_j_not_i = list(set(variabless_simple[indices[j]])-set(variabless_simple[indices[i]]))
#             vars_i_not_j.sort()
#             vars_j_not_i.sort()
#             print(i,j,len(variabless_simple[indices[i]]),len(variabless_simple[indices[j]]),len(constantss[indices[i]]),len(constantss[indices[j]]),sim_overlap,sim_jaccard,vars_i_not_j,vars_j_not_i)
#     print()

models_to_keep = ['T-Cell Signaling 2006_16464248.txt',
                  '27765040_tabular.txt',
                  'ErbB (1-4) Receptor Signaling_23637902.txt',
                  'HCC1954 Breast Cell Line Long-term ErbB Network_24970389.txt',
                  'T-LGL Survival Network 2011_22102804.txt',
                  'Predicting Variabilities in Cardiac Gene_26207376.txt',
                  'Lymphopoiesis Regulatory Network_26408858.txt',
                  'Lac Operon_21563979.txt',
                  'MAPK Cancer Cell Fate Network_24250280.txt',
                  'Septation Initiation Network_26244885.txt',
                  '29632237.txt',
                  '25063553_OR_OR.txt',
                  '19622164_TGF_beta1.txt',
                  '23658556_model_10.txt',
                  '23169817_high_dna_damage.txt',
                  '28426669_ARF10_greater_ARF5.txt',
                  '21450717_model_5_2.txt',
                  'Guard Cell Abscisic Acid Signaling_16968132.txt',
                  'FGF pathway of Drosophila Signaling Pathways_23868318.txt',
                  'Death Receptor Signaling_20221256.txt'
                  ]

models_to_exclude_manually_because_similar_from_same_PID = ['Trichostrongylus retortaeformis_22253585.txt',
                                                            'Bordetella bronchiseptica_22253585.txt']


def load_models_included_in_meta_analysis(max_degree=12,max_N=1000,similarity_threshold=0.9,folders=['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/'],models_to_keep=[],models_to_exclude_manually_because_similar_from_same_PID=[]):
## load the database, choose low max_n for quick results and to only look at small models
    [Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded] = db.load_database(folders,max_degree=max_degree,max_N=max_N)
    #similar_sets_jaccard = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=similarity_threshold,USE_JACCARD = True,models_to_keep=models_to_keep)[-1]
    #similar_sets_overlap = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=0.9,USE_JACCARD = False,models_to_keep=models_to_keep)[-1]
    #Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_excluded,similar_sets_jaccard = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=similarity_threshold,USE_JACCARD = True,models_to_keep=models_to_keep)
    Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_excluded,similar_sets = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=similarity_threshold,USE_JACCARD = False,models_to_keep=models_to_keep,models_to_exclude_manually_because_similar_from_same_PID=models_to_exclude_manually_because_similar_from_same_PID)
    n_variables = np.array(list(map(len,variabless)))
    n_constants = np.array(list(map(len,constantss)))
    return Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_excluded,models_not_loaded,similar_sets,n_variables,n_constants,max_degree

def nice_excel_table(data,xlabel,ylabel,cmap,upto,title):
    import xlsxwriter
    def rgba_to_hex(rgba):
        r,g,b,_ = list(map(int,np.array(list(rgba))*255))
        return '#%02x%02x%02x' % (r,g,b)
    
    workbook = xlsxwriter.Workbook(title)
    worksheet = workbook.add_worksheet()
    
    wrapallbutbottom = workbook.add_format({'text_wrap': True, 'top':True, 'bottom':False, 'left':True,'right':True, 'align': 'vcenter', 'valign': 'center', 'font_name':'Arial','font_size':10,'bold':1}) 
    wrapallbutright = workbook.add_format({'text_wrap': True, 'top':True, 'bottom':True, 'left':True,'right':False, 'align': 'vcenter', 'valign': 'center',  'rotation':90, 'font_name':'Arial','font_size':10,'bold':1}) 
    wrapallbuttop = workbook.add_format({'text_wrap': True, 'top':False, 'bottom':True, 'left':True,'right':True, 'align': 'vcenter', 'valign': 'center','font_name':'Arial','font_size':10,'bold':1}) 
    wrapallbutleft = workbook.add_format({'text_wrap': True, 'top':True, 'bottom':True, 'left':False,'right':True, 'align': 'vcenter', 'valign': 'center', 'font_name':'Arial','font_size':10,'bold':1}) 
    
    worksheet.merge_range(0,2,0,upto+2,xlabel,wrapallbutbottom)
    worksheet.merge_range(2,0,upto+2,0,ylabel,wrapallbutright)
    for i in range(upto+1):
        worksheet.write(1,2+i,(i),wrapallbuttop)    
        worksheet.write(2+i,1,(i),wrapallbutleft)
    
    for i in range(upto+1):
        for j in range(upto+1): 
            proportion = data[i,j]/sum(data[i,:])
            bg_color = rgba_to_hex(cmap(proportion)) if proportion>0 else 'FFFFFF'
            worksheet.write(2+i,2+j,data[i,j] if i>=j else '',workbook.add_format({'bg_color':bg_color,'font_color': 'black' if proportion<0.75 else 'white','top':i==0,'bottom':i==upto,'left':j==0,'right':j==upto,'font_name':'Arial','font_size':8,'align': 'vcenter', 'valign': 'center','bold':1}))
    worksheet.set_column(0,2+upto,4)
    worksheet.set_margins(top=0.1,left=0.1,right=0.1,bottom=0.1)
    workbook.close()

## General analyses

#interesting observations
#1) degree vs essential degree
def print_basic_summary_and_get_average_degree_per_model(Fs,Is,degrees,degrees_essential,variabless,constantss,N,n_variables,n_constants,max_degree):
    def summary(vec):
        return ', '.join(list(map(str,[np.min(vec),np.mean(vec),np.median(vec),np.max(vec)])))
    total_genes = sum(n_variables)
    vec_mean = []
    vec_total = []
    vec_max = []
    avg_degrees = []
    avg_essential_degrees = []
    for i in range(N):
        a=np.array(degrees[i][:n_variables[i]])
        b=np.array(degrees_essential[i][:n_variables[i]])
        dummy = np.array([np.array(el1)-np.array(el2) for el1,el2 in zip(a[a<=max_degree],b[a<=max_degree])])
        vec_mean.append(np.mean(dummy>0))
        vec_total.append(np.sum(dummy>0))
        vec_max.append(np.max(dummy))
        avg_degrees.append(np.mean(a))
        avg_essential_degrees.append(np.mean(b))
    

    #basic observations
    print( 'total genes:',total_genes )
    print( 'total constants:',sum(map(len,constantss)) )
    print( 'number of models with constants:', sum(np.array(list(map(len,constantss)))>0))
    print( 'number of genes per model (min, mean, median, max):',summary(n_variables) )
    print( 'average degree of models (min, mean, median, max):',summary(avg_degrees) )
    print( 'average essential degree of models (min, mean, median, max):',summary(avg_essential_degrees) )
    
    print( 'number of models that contain non-essential genes:', sum([el>0 for el in vec_total]))
    print('proportion of genes regulated by non-essential genes:',sum(vec_total),'/',total_genes,': ',sum(vec_total)/total_genes)

    return avg_degrees,avg_essential_degrees

def plot_degree_vs_actual_degree(degrees,degrees_essential,N,n_variables,max_degree,figsize=(3.43,3.1)):
    f,ax=plt.subplots(figsize=figsize)
    for i in range(N):
        a=np.array(degrees[i][:n_variables[i]])
        b=np.array(degrees_essential[i][:n_variables[i]])
        ax.plot(a[a<=max_degree],b[a<=max_degree],'ko',alpha=0.5)
    ax.set_xlabel('Number of inputs')
    ax.set_ylabel('Number of essential inputs')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim([ax.get_ylim()[0],ax.get_xlim()[1]])
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.savefig('degree_vs_actual_degree_N%i.pdf' % N,bbox_inches = "tight")    

def plot_size_vs_avg_connectivity(n_variables,avg_essential_degrees,N,figsize=(3.43,3.1)):
    import scipy.stats as stats
    x=np.array(n_variables)
    y = np.array(avg_essential_degrees)
    m,b,r,p,stderr = stats.linregress(np.log(x), y)
    r,p = stats.spearmanr(x, y)
    f,ax=plt.subplots(figsize=figsize)
    ax.semilogx(n_variables,avg_essential_degrees,'ko',markersize=8,alpha=0.5)
    ax.set_xlabel('Number of genes')
    ax.set_ylabel('Average essential degree')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.semilogx(x, b + m * np.log(x), 'x')
    (x1,x2) = ax.get_xlim()
    (y1,y2) = ax.get_ylim()
    #ax.semilogx([x1,x2],b+m*np.log([x1,x2]),'--',color='r',lw=2)
    ax.set_xlim([x1,x2])
    ax.text(np.exp(np.log(x1)+0.7*(np.log(x2)-np.log(x1))),y1+0.87*(y2-y1),(r'r$=%s$' % str(np.round(r,2))) +'\n'+('p$=%s$' % str(np.round(p,2))),color='red',va='center',ha='left')
    #print(stats.pearsonr(x,y))
    plt.gcf().subplots_adjust(bottom=0.2,left=0.2)
    plt.savefig('size_vs_avg_essential_degrees_N%i.pdf' % N,bbox_inches = "tight")

def plot_prop_pos_vs_avg_connectivity(prop_pos_separate,avg_essential_degrees,N,figsize=(3.43,3.1)):
    import scipy.stats as stats
    y=np.array(prop_pos_separate)
    x = np.array(avg_essential_degrees)
    m,b,r,p,stderr = stats.linregress(np.log(x), y)
    r,p = stats.spearmanr(x, y)
    f,ax=plt.subplots(figsize=figsize)
    ax.plot(x,y,'ko',markersize=8,alpha=0.5)
    ax.set_ylabel('Proportion activating edges')
    ax.set_xlabel('Average essential degree')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.semilogx(x, b + m * np.log(x), 'x')
    (x1,x2) = ax.get_xlim()
    (y1,y2) = ax.get_ylim()
    #ax.semilogx([x1,x2],b+m*np.log([x1,x2]),'--',color='r',lw=2)
    ax.set_xlim([x1,x2])
    ax.set_ylim([0,1.03])
    ax.text(np.exp(np.log(x1)+0.84*(np.log(x2)-np.log(x1))),y1+0.74*(y2-y1),(r'r$=%s$' % str(np.round(r,2))) +'\n'+('p$=%s$' % str(np.round(p,5))),color='red',va='center',ha='left')
    #print(stats.pearsonr(x,y))
    plt.gcf().subplots_adjust(bottom=0.2,left=0.2)
    plt.savefig('avg_essential_degrees_vs_prop_pos_N%i.pdf' % N,bbox_inches = "tight")


def plot_number_of_genes_vs_constants(n_variables,n_constants,N,figsize=(3.43,3.1)):
    import scipy.stats as stats
    f,ax=plt.subplots(figsize=figsize)
    dummy_n_constants = np.array(n_constants,dtype=float)
    dummy_n_constants[n_constants==0] = 0.5
    ax.loglog(n_variables,dummy_n_constants,'ko',markersize=8,alpha=0.5)
    #m,b,r,p,stderr = stats.linregress(n_variables, n_constants)
    (x1,x2) = ax.get_xlim()
    (y1,y2) = ax.get_ylim()
    #ax.semilogx([x1,x2],b+m*np.array([x1,x2]),'--',color='r',lw=2)
    r,p = stats.spearmanr(n_variables,n_constants)
    ax.set_xlim([x1,x2])
    ax.text(np.exp(np.log(x1)+0.03*(np.log(x2)-np.log(x1))),np.exp(np.log(y1)+0.87*(np.log(y2)-np.log(y1))),(r'r$=%s$' % str(np.round(r,2))) +'\n'+('p$=%s$' % str(np.round(p,2) if p>0.01 else '{:.0e}'.format(p))),color='red',va='center',ha='left')
    ax.set_xlabel('Number of genes')
    ax.set_ylabel('Number of external parameters')
    ax.set_yticks([0.5,1,10,100])
    ax.set_yticklabels(list(map(str,[0,1,10,100])))
    ax.spines['left'].set_visible(False)
    ax.plot([x1,x1],[y1,0.5],'k',lw=1)
    ax.plot([x1,x1],[1,y2],'k',lw=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)    
    minoryticks = np.array(ax.get_yticks(minor=True))
    minoryticks = minoryticks[minoryticks>1]
    minoryticks = minoryticks[minoryticks<y2]
    ax.set_yticks(minoryticks,minor=True)
    #kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    #d = .02
    #ax.plot((-d, +d), (0.1-d, 0.1+d), **kwargs)
    ax.set_ylim([y1,y2]) 
    plt.gcf().subplots_adjust(bottom=0.2,left=0.2)
    plt.savefig('number_of_genes_vs_constants_nice_N%i.pdf' % N,bbox_inches = "tight")

def plot_summary_all_networks(n_variables,n_constants,N):
    sorted_n_variables = np.array(sorted(zip(n_variables,n_constants),key=lambda x: x[0])).T
    f,ax=plt.subplots(figsize=(7.01,2.21))
    ind = np.arange(N)
    ax.bar(ind,sorted_n_variables[0],width=0.7,color=[0.8,0,0],label='genes (%s total)' % str(sum(n_variables)))
    ax.bar(ind,sorted_n_variables[1],width=0.7,bottom=sorted_n_variables[0],color=[0.2,0.2,0.9],label='external parameters (%s total)' % str(sum(n_constants)))
    ax.set_xlim([-1,N])
    ax.set_ylabel('Number of nodes')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([])
    ax.set_xlabel(str(N)+' networks sorted by number of genes')
    ax.legend(loc=2,frameon=False)
    #ax.text(N/2.7,300,str(N)+' networks\n'+str(sum(n_variables))+' genes\n'+str(sum(n_constants))+' external parameters',verticalalignment='center',horizontalalignment='left')
    plt.savefig('summary_all_networks_N%i.pdf' % N,bbox_inches = "tight")

def print_type_of_function(Fs):
    ## proportion of increasing vs decreasing fcts
    res = []
    all_res = []
    for i in range(len(Fs)):
        F = Fs[i]
        res.append([])
        for f in F:
            if len(f)==0: #happens if actual degree_f > max_degree
                continue
            dummy = can.is_monotonic(f,True)[1]
            res[-1].append(dummy)
            all_res.extend(dummy)
    print(pd.Series(all_res).value_counts())

def print_common_variables_across_networks(variabless):
    ## Which variables occur in many networks?
    variabless_simple = [[el.lower().replace('_','').replace('.','').replace('kappa','k') for el in variables] for variables in variabless]
    
    
    all_variables = []
    for i in range(len(variabless)):
        #all_variables.extend(list(map(str.lower,variabless[i])))
        all_variables.extend(variabless_simple[i])
    all_variables = np.array(all_variables)
    print(pd.Series(all_variables).value_counts()[:50])
    #common_variables = list(pd.Series(all_variables).value_counts()[:20].index)

def plot_in_and_out_degree(Is,n_variables,degrees,N,figsize=(3.43,2.6)):
    ## In-degree and Out-degree distribution, also for common variables
    all_degrees = []
    for i in range(N):
        all_degrees.extend(degrees[i][:n_variables[i]])
    all_degrees = np.array(all_degrees)
    
    all_outdegrees = []
    outdegrees = []
    all_outdegrees_prop = []
    for i in range(N):
        outdegrees.append([0 for j in range(n_variables[i])])
        for regulators in Is[i][:n_variables[i]]:
            for regulator in regulators:
                try:
                    outdegrees[i][regulator] += 1
                except IndexError: #if regulator is a constant
                    pass
        all_outdegrees.extend(outdegrees[i])
        mean_outdegree = np.mean(outdegrees[i])
        all_outdegrees_prop.extend([el/mean_outdegree for el in outdegrees[i]])
    all_outdegrees = np.array(all_outdegrees)
    all_outdegrees_prop = np.array(all_outdegrees_prop)
    
    #outdegrees_common_variables = [np.mean(all_outdegrees[all_variables==common_variables[i]]) for i in range(20)]
    #outdegrees_prop_common_variables = [np.mean(all_outdegrees_prop[all_variables==common_variables[i]]) for i in range(20)]
    
    dummy = pd.Series(all_outdegrees).value_counts()
    all_outdegrees_index = np.array(dummy.index)
    all_outdegrees_values = np.array(dummy)
    # f,ax = plt.subplots()
    # for outdegree,count in zip(all_outdegrees_index,all_outdegrees_values):
    #     if outdegree>0:
    #         ax.loglog([outdegree],[count],'ko')
    # ax.set_xlabel('out-degree')
    # ax.set_ylabel('total count')
    # plt.savefig('distribution_outdegree_N%i.pdf' % len(Fs),bbox_inches = "tight")
    
    dummy = pd.Series(all_degrees).value_counts()
    all_degrees_index = np.array(dummy.index)
    all_degrees_values = np.array(dummy)
    # f,ax = plt.subplots()
    # for degree,count in zip(all_degrees_index,all_degrees_values):
    #     if degree>0:
    #         ax.loglog([degree],[count],'ko')
    # ax.set_xlabel('in-degree')
    # ax.set_ylabel('total count')
    # plt.savefig('distribution_indegree_N%i.pdf' % len(Fs),bbox_inches = "tight")
    
    # #in degree and out degree in one plot
    # f,ax = plt.subplots()
    # ax.loglog(all_degrees_index,all_degrees_values/sum(all_degrees_values),'ro',label='in-degree')
    # #ax.loglog(all_outdegrees_index[all_outdegrees_index>0],all_outdegrees_values[all_outdegrees_index>0]/sum(all_outdegrees_values[all_outdegrees_index>0]),'k*',label='out-degree')
    # ax.loglog(all_outdegrees_index[all_outdegrees_index>0],all_outdegrees_values[all_outdegrees_index>0]/sum(all_outdegrees_values),'k*',label='out-degree')
    # ax.set_xlabel('degree')
    # ax.set_ylabel('proportion')
    # ax.legend(loc='best')
    # plt.savefig('distribution_in_and_outdegree_N%i.pdf' % len(Fs),bbox_inches = "tight")
    
    #in degree and out degree in one plot, nice
    f,ax = plt.subplots(figsize=figsize)
    ax.loglog(all_degrees_index,all_degrees_values/sum(all_degrees_values),'ro',label='in-degree')
    all_outdegrees_index_dummy = np.array(all_outdegrees_index,dtype=float)
    all_outdegrees_index_dummy[all_outdegrees_index==0] = 0.5
    #ax.loglog(all_outdegrees_index[all_outdegrees_index>0],all_outdegrees_values[all_outdegrees_index>0]/sum(all_outdegrees_values[all_outdegrees_index>0]),'k*',label='out-degree')
    ax.loglog(all_outdegrees_index_dummy,all_outdegrees_values/sum(all_outdegrees_values),'k*',label='out-degree')
    (x1,x2) = ax.get_xlim()
    (y1,y2) = ax.get_ylim()
    y1=1e-4
    y2=1
    ax.set_xlabel('Degree')
    ax.set_ylabel('Proportion')
    ax.legend(loc='best',frameon=False)
    ax.set_xticks([0.5,1,2,5,10,20,50])
    ax.set_xticklabels(list(map(str,[0,1,2,5,10,20,50])))
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.plot([x1,0.5],[y1,y1],'k',lw=1,clip_on=False)
    ax.plot([1,x2],[y1,y1],'k',lw=1,clip_on=False)
    minorxticks = np.array(ax.get_xticks(minor=True))
    minorxticks = minorxticks[minorxticks>1]
    minorxticks = minorxticks[minorxticks<x2]
    ax.set_xticks(minorxticks,minor=True)
    ax.set_xlim([x1,x2]) 
    ax.set_ylim([y1,y2])
    #ax2 = ax.twinx()
    #ax2.loglog(all_outdegrees_index_dummy,all_outdegrees_values,'k*')
    #ax2.set_ylabel('total count')
    #ax2.spines['bottom'].set_visible(False)
    #ax2.xticks['bottom'].set_visible(False)
    #ax.set_ylim([1e-4,1])
    plt.savefig('distribution_in_and_outdegree_N%i_nice.pdf' % N,bbox_inches = "tight")

def compute_type_of_each_regulation(Fs):
    type_of_each_regulation = []
    for i,F in enumerate(Fs):
        for f in F:
            if len(f)==0: #happens if actual degree_f > max_degree
                continue
            (NONDEGENERATED,monotonic) = can.is_monotonic(f,True)
            #nrs = associate_number_monotonic(*monotonic)
            n_essential=sum([el!='not essential' for el in monotonic])
            n = len(monotonic)
            for el in monotonic:
                type_of_each_regulation.append([i,n_essential,n,el])
    type_of_each_regulation = pd.DataFrame(type_of_each_regulation,columns = ['ID','essential_degree','degree','type'])
    return type_of_each_regulation

def compute_type_of_each_regulation_excluding_constants(Fs,n_variables):
    type_of_each_regulation = []
    for i,F in enumerate(Fs):
        for f in F[:n_variables[i]]:
            if len(f)==0: #happens if actual degree_f > max_degree
                continue
            (NONDEGENERATED,monotonic) = can.is_monotonic(f,True)
            #nrs = associate_number_monotonic(*monotonic)
            n_essential=sum([el!='not essential' for el in monotonic])
            n = len(monotonic)
            for el in monotonic:
                type_of_each_regulation.append([i,n_essential,n,el])
    type_of_each_regulation = pd.DataFrame(type_of_each_regulation,columns = ['ID','essential_degree','degree','type'])
    return type_of_each_regulation


def compute_type_of_each_regulation_by_model(type_of_each_regulation,N):
    type_of_each_regulation_by_model = np.zeros((N,4),dtype=int)
    for i in range(N):
        for j,t in enumerate(['increasing','decreasing','not monotonic','not essential']):
            type_of_each_regulation_by_model[i,j] = sum(type_of_each_regulation.type[type_of_each_regulation.ID==i]==t)
    return type_of_each_regulation_by_model


def plot_activators_vs_inhibitors_vs_degree(type_of_each_regulation,N,figsize=(3.43,3.1),SHOW_N=False):
    color_unknown = 'black'#[0.5,0.5,0.5]
    color_neg = 'orange'#[0.7,0.7,1]
    color_pos = 'blue'#[1,0.7,0.7]
    
    types = ['increasing','decreasing','not monotonic']
    
    n_max = 7
    colors = [color_pos,color_neg,color_unknown]
    #f,ax = plt.subplots()
    x = np.arange(1,n_max+1)
    type_of_each_regulation_aggregated = []
    for i,value in enumerate(types):
        y = np.array([sum(type_of_each_regulation.type[type_of_each_regulation.essential_degree==n]==value) for n in range(1,n_max+1)])
        type_of_each_regulation_aggregated.append(y)
        #ax.semilogy(x[y>0],y[y>0] ,'x-',color=colors[i])
    #ax.legend(['increasing','decreasing','not monotonic'],loc='best')
    type_of_each_regulation_aggregated = np.array(type_of_each_regulation_aggregated)
    
    f,ax = plt.subplots(figsize=figsize)
    type_of_each_regulation_aggregated_prop = type_of_each_regulation_aggregated/sum(type_of_each_regulation_aggregated,0)
    for i,label in enumerate(types):
        ax.bar(x,type_of_each_regulation_aggregated_prop[i,:],bottom=np.sum(type_of_each_regulation_aggregated_prop[:i,:],0),color=colors[i],alpha=0.3)
    for i,val in enumerate(x):
        ax.text(val,1.05,'n='+str(np.sum(type_of_each_regulation_aggregated[:,i])),va='center',ha='center',clip_on=False,fontsize=7)
    #ax.legend(['activation','inhibition','conditional'],bbox_to_anchor=(1.05, 0.65))
    #ax.legend(['activation','inhibition','conditional'],loc=8,ncol=1)
    ax.legend(['activating','inhibitory','conditional'],loc=8,ncol=1,title='Type of regulation')
    ax.set_ylim([0,1])
    ax.set_xticks(list(range(1,n_max+1)))
    ax.set_xlabel('Number of essential inputs')
    ax.set_ylabel('Proportion')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig('activators_vs_inhibitors_vs_degree_N%s.pdf' % str(N),bbox_inches = "tight")












## Canalization

def get_canalizing_depths(Fs,degrees_essential,n_variables,N):
    ## calculate the "degree" of canalization for each function
    depths,layer_structures,all_ncfs,layer_structures_ncf = [],[],[],[]
    ks_per_n = np.zeros((max_degree+15,max_degree+15))
    for i in range(N):
        F = Fs[i]
        depths.append([])
        layer_structures.append([])
        layer_structures_ncf.append([])
        for j in range(n_variables[i]):
            f = F[j]
            if len(f)==0: #happens if actual degree_f > max_degree
                depth = np.nan
                layer_structure = []
            else:
                (n_f,depth,can_inputs,can_outputs,corefunction) = can.get_canalizing_depth_inputs_outputs_corefunction(f)
                layer_structure = can.get_layerstructure_given_canalizing_outputs_and_corefunction(can_outputs,corefunction,n_f)
            depths[-1].append(depth)
            layer_structures[-1].append(layer_structure)
            if n_f == depth:
                layer_structures_ncf[-1].append(', '.join(list(map(str,layer_structure))))
        all_ncfs.append(sum(depths[-1]) == sum(degrees_essential[i][:n_variables[i]]))
        for k,n in zip(depths[-1],degrees_essential[i]):
            if not np.isnan(k):
                ks_per_n[n,k] += 1
    return depths,layer_structures,all_ncfs,layer_structures_ncf,ks_per_n
            
def write_excel_files_canalizing_depth(ks_per_n,N):
    
    ## various observations about canalization
    ks_per_n_prop = np.copy(ks_per_n)
    for i in range(len(ks_per_n_prop)):
        sum_row = sum(ks_per_n_prop[i,:])
        if sum_row>0:
            ks_per_n_prop[i,:] = ks_per_n_prop[i,:] *100./sum_row
    
    ks_per_n_prop_greaterthan = np.copy(ks_per_n_prop)
    for i in range(len(ks_per_n_prop)):
        ks_per_n_prop_greaterthan[i,:] = np.append(np.cumsum(ks_per_n_prop[i,:0:-1]),ks_per_n_prop[i,0])[::-1]
    
    pd.DataFrame(ks_per_n,index=['n=%i' % n for n in range(len(ks_per_n))],columns=['k=%i' % n for n in range(len(ks_per_n))]).to_excel('ks_per_n_%inetworks.xlsx' % N)
    pd.DataFrame(ks_per_n_prop,index=['n=%i' % n for n in range(len(ks_per_n))],columns=['k=%i' % n for n in range(len(ks_per_n))]).to_excel('ks_per_n_proportions_%inetworks.xlsx' % N)
    pd.DataFrame(ks_per_n_prop_greaterthan,index=['n=%i' % n for n in range(len(ks_per_n))],columns=['k=0']+['k>=%i' % n for n in range(1,len(ks_per_n))]).to_excel('ks_per_n_proportions_greaterthan_%inetworks.xlsx' % N)

    from matplotlib import cm
    nice_excel_table(ks_per_n,xlabel='canalizing depth',ylabel='number of essential inputs',cmap=cm.Blues,upto=20,title='ks_per_n_%inetworks_mod.xlsx' % N)


    upto = 10
    nsim = 1000
    ks_per_n_random = np.zeros((upto+1,upto+1))
    for n in range(upto+1):
        for _ in range(nsim):
            F = can.random_non_degenerated_function(n)
            k = can.get_canalizing_depth_inputs_outputs_corefunction(F)[1]
            ks_per_n_random[n,k] += 1
                
    ks_per_n_random_prop = np.copy(ks_per_n_random)
    for i in range(len(ks_per_n_random_prop)):
        sum_row = sum(ks_per_n_random_prop[i,:])
        if sum_row>0:
            ks_per_n_random_prop[i,:] = ks_per_n_random_prop[i,:] *100./sum_row
    
    ks_per_n_random_prop_greaterthan = np.copy(ks_per_n_random_prop)
    for i in range(len(ks_per_n_random_prop)):
        ks_per_n_random_prop_greaterthan[i,:] = np.append(np.cumsum(ks_per_n_random_prop[i,:0:-1]),ks_per_n_random_prop[i,0])[::-1]
    
    pd.DataFrame(ks_per_n_random,index=['n=%i' % n for n in range(len(ks_per_n_random))],columns=['k=%i' % n for n in range(len(ks_per_n_random))]).to_excel('ks_per_n_random_nsim%i.xlsx' % nsim)
    pd.DataFrame(ks_per_n_random_prop,index=['n=%i' % n for n in range(len(ks_per_n_random))],columns=['k=%i' % n for n in range(len(ks_per_n_random))]).to_excel('ks_per_n_proportions_random_nsim%i.xlsx' % nsim)
    pd.DataFrame(ks_per_n_random_prop_greaterthan,index=['n=%i' % n for n in range(len(ks_per_n_random))],columns=['k=0']+['k>=%i' % n for n in range(1,len(ks_per_n_random))]).to_excel('ks_per_n_proportions_greaterthan_random_nsim%i.xlsx' % nsim)
    
    nice_excel_table(ks_per_n_random,xlabel='canalizing depth',ylabel='number of essential inputs',cmap=cm.Reds,upto=upto,title='ks_per_n_random_%inetworks_mod.xlsx' % N)




def get_number_of_layers(Fs):
    res = []
    for ii,F in enumerate(Fs):
        for f in F:
            if len(f)>0: #if degree(f) > max_degree, f==[]
                (n,k,can_inputs,can_outputs,corefunction) = can.get_canalizing_depth_inputs_outputs_corefunction(f)
                if k>-1:
                    kis = can.get_layerstructure_given_canalizing_outputs_and_corefunction(can_outputs,corefunction,n)
                    r = len(kis)
                else:
                    kis = []
                    r = 0
                res.append([ii,n,k,can_inputs,can_outputs,kis,r,np.mean(f)])
        
    res = np.array(res)
    
    upto = 10
    obs_counts = np.zeros((upto-3+1,upto))
    exp_counts = np.zeros((upto-3+1,upto))
    for nk in range(3,upto+1):
        a = pd.Series(res[np.bitwise_and(res[:,1]==nk,res[:,2]==nk)][:,5])
        print(a.value_counts())
        
        a = pd.Series(res[np.bitwise_and(res[:,1]==nk,res[:,2]==nk)][:,6])
        print(a.value_counts())
        dummy = a.value_counts()
        for key in dummy.keys():
            obs_counts[nk-3,key-1]  = dummy[key]
        
        nsim = 1000
        res_rnd = []
        for i in range(nsim):
            res_rnd.append(can.get_layerstructure_given_canalizing_outputs_and_corefunction(can.get_canalizing_depth_inputs_outputs_corefunction(can.random_k_canalizing(nk,nk,True))[3],[0],nk))
        a = pd.Series(res_rnd)
        print(a.value_counts())
        
        a = pd.Series(list(map(len,res_rnd)))
        print(a.value_counts())
        dummy = a.value_counts()
        for key in dummy.keys():
            exp_counts[nk-3,key-1]  = dummy[key]
            
    pd.DataFrame(obs_counts,index=['n = '+el for el in list(map(str,range(3,upto+1)))],columns=list(map(str,range(1,upto+1)))).to_excel('observed_number_of_layers_n%i.xlsx' % len(Fs))
    pd.DataFrame(exp_counts,index=['n = '+el for el in list(map(str,range(3,upto+1)))],columns=list(map(str,range(1,upto+1)))).to_excel('expected_number_of_layers_nsim%i.xlsx' % (nsim))
    return res

def get_canalizing_input_and_output_distribution_of_NCFs(res_from_get_number_of_layers,N,degree=3):
    res = res_from_get_number_of_layers
    
    nk = degree
    a = pd.Series([list(el[:-1]) for el in res[np.bitwise_and(res[:,1]==nk,res[:,2]==nk)][:,4]])
    ex = a.value_counts()
    pd.DataFrame(np.c_[list(ex.index),ex.values],columns = ['b'+str(i+1) for i in range(nk-1)]+['count']).to_excel('can_output_distribution_for_ncfs_with_n%i.xlsx' % nk)
    
    nk = degree
    a = pd.Series([list(el[:-1]) for el in res[np.bitwise_and(res[:,1]==nk,res[:,2]==nk)][:,3]])
    ex = a.value_counts()
    pd.DataFrame(np.c_[list(ex.index),ex.values],columns = ['a'+str(i+1) for i in range(nk-1)]+['count']).to_excel('can_input_distribution_for_ncfs_with_n%i.xlsx' % nk)
    
    f,ax = plt.subplots()
    for nk in range(1,11):
        ax.bar([nk],[np.mean(res[np.bitwise_and(res[:,1]==nk,res[:,2]==nk)][:,7])],color='k')
    ax.set_xlabel('Degree')
    ax.set_ylabel('Average proportion of 1s')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.savefig('avg1s_vs_n_N%i.pdf' % N ,bbox_inches = "tight")
    
    
    

#layer structures
def get_count_of_specific_NCF_layer_structures(layer_structures_ncf,N):
    layer_structures_ncf_flat = pd.Series([entry for el in layer_structures_ncf for entry in el])
    values = layer_structures_ncf_flat.value_counts()
    index = values.index
    values = values.values
    ns = np.array([sum(map(int,el.split(', '))) if len(el)>0 else 0 for el in index])
    
    f = open('layers_structure_counts_N%i.txt' % N,'w')
    for n in range(2,7):
        f.write('n=%i\n' % n)
        for a,b in zip(index[ns==n],values[ns==n]):
            f.write('%i\t%s\n' % (b,a))
        f.write('\n\n')
    f.close()


# for n in range(2,7):
#     res = []
#     for i in range(1000):
#         f = can.random_k_canalizing(n,n)
#         (n,k,inputs,outputs,corefunction) = can.get_canalizing_depth_inputs_outputs_corefunction(f)
#         layers = can.get_layer_structure_given_outputs_corefunction(outputs,corefunction,n)
#         res.append(', '.join(list(map(str,layers))))
#     res = pd.Series(res)
#     values = res.value_counts()
#     index = values.index
#     values = values.values
#     ns = np.array([sum(map(int,el.split(', '))) for el in index])

#     for a,b in zip(index[ns==n],values[ns==n]):
#         print()
#     print()



def get_symmetry_groups_of_networks(Fs,degrees,degrees_essential,N,upto=11):
    import itertools
    #calculate symmetry groups
    nr_symmetry_groups =[]
    ss_per_n = np.zeros((upto+1,upto+1))
    for i in range(N):
        F = Fs[i]
        nr_symmetry_groups.append([])
        for j,f in enumerate(F):
            if len(f)==0 or degrees_essential[i][j]>upto: #happens if actual degree_f > max_degree
                number = np.nan
            else:
                if degrees[i][j]==degrees_essential[i][j]:
                    number = len(can.get_symmetry_groups(f))
                else:
                    indices_essential_variables = can.get_essential_variables(f)
                    indices_non_essential_variables = list(set(range(degrees[i][j])) - set(indices_essential_variables))
                    bool_list = np.array(list(itertools.product([0, 1], repeat=degrees[i][j])))
                    f_ess = np.array(f)[np.array(db.find_all_indices(np.sum(bool_list[:,np.array(indices_non_essential_variables)],1)==0,True))]
                    number = len(can.get_symmetry_groups(f_ess))
                ss_per_n[ degrees_essential[i][j] , number ] += 1
            nr_symmetry_groups[-1].append(number)
    pd.DataFrame(ss_per_n,index=['n=%i' % n for n in range(len(ss_per_n))],columns=['s=%i' % n for n in range(len(ss_per_n))]).to_excel('nr_symmetry_groups_per_n_%inetworks.xlsx' % N)
    from matplotlib import cm
    nice_excel_table(ss_per_n,xlabel='number of symmetry groups',ylabel='number of essential inputs',cmap=cm.Blues,upto=upto,title='ss_per_n_%inetworks_mod.xlsx' % N)
    return ss_per_n

def get_random_expectations_symmetry_groups(ks_per_n,upto=11,nsim = 1000):
    import itertools
    ss_per_n_random = np.zeros((upto+1,upto+1))
    for n in range(1,upto+1):
        bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
        for _ in range(nsim):
            F = can.random_non_degenerated_function(n)
            s = len(can.get_symmetry_groups(F,bool_list))
            ss_per_n_random[n,s] += 1
    pd.DataFrame(ss_per_n_random,index=['n=%i' % n for n in range(len(ss_per_n_random))],columns=['s=%i' % n for n in range(len(ss_per_n_random))]).to_excel('nr_symmetry_groups_per_n_random_nsim%i.xlsx' % (nsim))
    
    ss_per_n_random_canalizing = np.zeros((upto+1,upto+1))
    for n in range(1,upto+1):
        bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
        for _ in range(nsim):
            F = can.random_k_canalizing(n,1,False,bool_list)
            s = len(can.get_symmetry_groups(F,bool_list))
            ss_per_n_random_canalizing[n,s] += 1
    pd.DataFrame(ss_per_n_random_canalizing,index=['n=%i' % n for n in range(len(ss_per_n_random_canalizing))],columns=['s=%i' % n for n in range(len(ss_per_n_random_canalizing))]).to_excel('nr_symmetry_groups_per_n_random_canalizing_nsim%i.xlsx' % (nsim))
    
    ss_per_n_random_NCF = np.zeros((upto+1,upto+1))
    for n in range(1,upto+1):
        bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
        for _ in range(nsim):
            F = can.random_k_canalizing(n,n,False,bool_list)
            s = len(can.get_symmetry_groups(F,bool_list))
            ss_per_n_random_NCF[n,s] += 1
    pd.DataFrame(ss_per_n_random_NCF,index=['n=%i' % n for n in range(len(ss_per_n_random_NCF))],columns=['s=%i' % n for n in range(len(ss_per_n_random_NCF))]).to_excel('nr_symmetry_groups_per_n_random_NCF_nsim%i.xlsx' % (nsim))
    
    ss_per_n_imputed = np.zeros((upto+1,upto+1))
    for n in range(1,upto+1):
        bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
        for k in np.random.choice(list(range(n+1)),nsim,replace=True,p=ks_per_n[n][:(n+1)]/sum(ks_per_n[n][:(n+1)])):
            F = can.random_k_canalizing(n,k,True,bool_list)
            s = len(can.get_symmetry_groups(F,bool_list))
            ss_per_n_imputed[n,s] += 1
    pd.DataFrame(ss_per_n_imputed,index=['n=%i' % n for n in range(len(ss_per_n_imputed))],columns=['s=%i' % n for n in range(len(ss_per_n_imputed))]).to_excel('nr_symmetry_groups_per_n_random_imputed_nsim%i.xlsx' % (nsim))

    from matplotlib import cm
    nice_excel_table(ss_per_n_random,xlabel='number of symmetry groups',ylabel='number of essential inputs',cmap=cm.Reds,upto=upto,title='ss_per_n_random_%inetworks_mod.xlsx' % N)
    nice_excel_table(ss_per_n_random_canalizing,xlabel='number of symmetry groups',ylabel='number of essential inputs',cmap=cm.Reds,upto=upto,title='ss_per_n_random_canalizing_%inetworks_mod.xlsx' % N)
    nice_excel_table(ss_per_n_random_NCF,xlabel='number of symmetry groups',ylabel='number of essential inputs',cmap=cm.Reds,upto=upto,title='ss_per_n_random_NCF_%inetworks_mod.xlsx' % N)
    nice_excel_table(ss_per_n_imputed,xlabel='number of symmetry groups',ylabel='number of essential inputs',cmap=cm.Reds,upto=upto,title='ss_per_n_imputed_%inetworks_mod.xlsx' % N)
    return ss_per_n_random,ss_per_n_random_canalizing,ss_per_n_random_NCF,ss_per_n_imputed

def plot_symmetry_distributions(ss_per_n,ss_per_n_imputed,ss_per_n_random,N):
    f,ax = plt.subplots(figsize=(4,3.1))
    from matplotlib import cm
    cmap = cm.tab10
    width=0.55
    height=1
    for j,data in enumerate([ss_per_n,ss_per_n_imputed,ss_per_n_random]):
        for i in range(2,ss_per_n.shape[0]):
            norm_data = data[i,:]/sum(data[i,:])
            ax.barh(np.arange(len(norm_data)),width=width*norm_data,left=4*i+j,color=cmap(j),height=height)
            ax.barh(np.arange(len(norm_data)),width=width*norm_data,left=4*i+j-width*norm_data,color=cmap(j),height=height)

    for j,(data,label) in enumerate(zip([ss_per_n,ss_per_n_imputed,ss_per_n_random],['observed','random,\ncanalizing depth matched','random'])):
        ax.barh(np.arange(len(norm_data)),width=width*norm_data,left=-100,color=cmap(j),height=height,label=label)
    ax.set_xticks(8+1+4*np.arange(ss_per_n.shape[0]-2))
    ax.set_xticklabels(list(map(str,range(2,ss_per_n.shape[0]))))
    ax.set_xlabel('Number of essential inputs')
    ax.set_ylabel('Number of symmetry groups')
    ax.set_xlim([8-1,4*10+3])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim([0.4,10.6])
    ax.legend(loc='upper left',frameon=False)
    plt.savefig('symmetry_groups%i.pdf' % (N),bbox_inches = "tight")
    






def get_canalizing_strengths(Fs):
    ## Canalizing strength - needs work!  
    can_strengths = []
    for F in Fs:
        can_strengths.append([])
        for f in F:
            len_f = len(f)
            if len_f>2 and len_f<33:
                can_strengths[-1].append(can.get_canalizing_strength(f))
    return can_strengths

def plot_canalizing_strength_of_noncanalizing_functions(N,max_degree,nsim=100,figsize=(3.43,3.1)):
    import itertools
    import matplotlib.patches as mpatches
    #canalizing strength of non-canalizing functions (n>=3)
    a,b = 3,6
    k = 0
    can_strengths = [[] for _ in range(max_degree)]
    bool_lists = [np.array(list(itertools.product([0, 1], repeat=n))) for n in range(1,b+1)]
    non_canalizing_functions = [[] for _ in range(max_degree)]
    for F in Fs:
        for f in F:
            if len(f)==0:
                continue
            (n_f,depth,can_inputs,can_outputs,corefunction) = can.get_canalizing_depth_inputs_outputs_corefunction(f)
            #n_f = can.nr_essential_variables(f)
            if depth == k and n_f>=a and n_f<=b:
                can_strengths[n_f-1].append( can.get_canalizing_strength(f,bool_lists[n_f-1])[0] )
                non_canalizing_functions[n_f-1].append(f)
    
    can_strengths_random = [[] for _ in range(max_degree)]
    non_canalizing_functions_random = [[] for _ in range(max_degree)]
    for n in range(a,b+1):
        for ii in range(nsim):
            f_random = can.random_k_canalizing(n,k,True)
            can_strengths_random[n-1].append( can.get_canalizing_strength(f_random,bool_lists[n-1])[0] )
            non_canalizing_functions_random[n-1].append(f_random)
    
    f,ax = plt.subplots(figsize=figsize)
    
    labels = []
    def add_label(violin, label):
        color = violin["bodies"][0].get_facecolor().flatten()
        labels.append((mpatches.Patch(color=color), label))
        
    positions = np.arange(0.5-0.42,(b-a+1)*2,2)
    data = can_strengths[a-1:b]
    add_label(ax.violinplot(data,positions=positions,showextrema=False,showmeans=True,widths=0.8), "Observed")    
    
    positions = np.arange(0.5+0.42,(b-a+1)*2,2)
    data = can_strengths_random[a-1:b]
    add_label(ax.violinplot(data,positions=positions,showextrema=False,showmeans=True,widths=0.8), "Random")    
    
    ax.set_xticks(np.arange(0.5,(b-a+1)*2,2))
    ax.set_xticklabels([str(degree) for degree in range(a,b+1)])
    ax.set_xlabel('Number of essential inputs')
    ax.set_ylabel('Canalizing strength')
    ax.set_ylim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(*zip(*labels),loc='best',frameon=False)
    plt.gcf().subplots_adjust(bottom=0.17,left=0.2)
    ax.set_title('')
    plt.savefig('canalizing_strength_nmin%i_nmax%i_k%i_N%i_nsim%i.pdf' % (a,b,k,N,nsim),bbox_inches = "tight")

def plot_effectiveness_of_noncanalizing_functions(N,max_degree,nsim=100,figsize=(3.43,3.1)):
    import matplotlib.patches as mpatches
    #canalizing strength of non-canalizing functions (n>=3)
    a,b = 3,6
    k = 0
    effectivenesses = [[] for _ in range(max_degree)]
    non_canalizing_functions = [[] for _ in range(max_degree)]
    for F in Fs:
        for f in F:
            if len(f)==0:
                continue
            (n_f,depth,can_inputs,can_outputs,corefunction) = can.get_canalizing_depth_inputs_outputs_corefunction(f)
            #n_f = can.nr_essential_variables(f)
            if depth == k and n_f>=a and n_f<=b:
                effectivenesses[n_f-1].append( can.get_input_redundancy(f) )
                non_canalizing_functions[n_f-1].append(f)
    
    effectivenesses_random = [[] for _ in range(max_degree)]
    non_canalizing_functions_random = [[] for _ in range(max_degree)]
    for n in range(a,b+1):
        for ii in range(nsim):
            f_random = can.random_k_canalizing(n,k,True)
            effectivenesses_random[n-1].append( can.get_input_redundancy(f_random) )
            non_canalizing_functions_random[n-1].append(f_random)
    
    f,ax = plt.subplots(figsize=figsize)
    
    labels = []
    def add_label(violin, label):
        color = violin["bodies"][0].get_facecolor().flatten()
        labels.append((mpatches.Patch(color=color), label))
        
    positions = np.arange(0.5-0.42,(b-a+1)*2,2)
    data = effectivenesses[a-1:b]
    add_label(ax.violinplot(data,positions=positions,showextrema=False,showmeans=True,widths=0.8), "Observed")    
    
    positions = np.arange(0.5+0.42,(b-a+1)*2,2)
    data = effectivenesses_random[a-1:b]
    add_label(ax.violinplot(data,positions=positions,showextrema=False,showmeans=True,widths=0.8), "Random")    
    
    ax.set_xticks(np.arange(0.5,(b-a+1)*2,2))
    ax.set_xticklabels([str(degree) for degree in range(a,b+1)])
    ax.set_xlabel('Number of essential inputs')
    ax.set_ylabel('Normalized input redundancy')
    ax.set_ylim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(*zip(*labels),loc='best',frameon=False)
    plt.gcf().subplots_adjust(bottom=0.17,left=0.2)
    ax.set_title('')
    plt.savefig('effectiveness_nmin%i_nmax%i_k%i_N%i_nsim%i.pdf' % (a,b,k,N,nsim),bbox_inches = "tight")


def compute_effectiveness_of_each_regulation(Fs,max_degree):
    try:
        import cana
    except ModuleNotFoundError:
        print("ERROR: module 'cana' not installed")
        return 
    effectiveness_of_each_regulation = []
    for i,F in enumerate(Fs):
        for f in F:
            if len(f) == 0 or np.log2(len(f))>max_degree: #happens if actual degree_f > max_degree
                continue
            n = int(np.log2(len(f)))
            edge_effectiveness = np.array(cana.boolean_node.BooleanNode(k=n,inputs=range(n),outputs=f).edge_effectiveness())
            essential_degree = sum(edge_effectiveness>0)
            for el in edge_effectiveness:
                effectiveness_of_each_regulation.append([i,essential_degree,n,el])
    effectiveness_of_each_regulation = pd.DataFrame(effectiveness_of_each_regulation,columns = ['ID','essential_degree','degree','effectiveness'])
    return effectiveness_of_each_regulation







## Dynamics - Derrida values

def get_derrida_values(Fs,Is,degrees,max_degree,nsim=1000):
    all_max_degrees = list(map(max,degrees))
    ## Derrida plot - some basic plots
    derridas = []
    for i in range(len(Fs)):
        if all_max_degrees[i]>max_degree:
            derridas.append(np.nan)
            continue
        F = Fs[i]
        I = Is[i]
        len_F = len(F)
        derrida = can.derrida_value(F,I,len_F,1,nsim=nsim)
        derridas.append(derrida)
    return derridas

def plot_derrida_values_vs_network_size(derridas,n_variables,max_degree,N):
    import scipy.stats as stats
    all_max_degrees = list(map(max,degrees)) 

    n_variables_restr = [n_variables[i] for i in range(N) if all_max_degrees[i]<=max_degree]
    derridas_restr = [derridas[i] for i in range(N) if all_max_degrees[i]<=max_degree]

    f,ax = plt.subplots()
    ax.semilogx(n_variables_restr,derridas_restr,'ko',alpha=0.5)
    ax.set_xlabel('Number of genes')
    ax.set_ylabel('Derrida value for a\n single perturbation')
    (x1,x2) = ax.get_xlim()
    (y1,y2) = ax.get_ylim()
    r,p = stats.spearmanr(n_variables_restr, derridas_restr)
    ax.text(np.exp(np.log(x1)+0.7*(np.log(x2)-np.log(x1))),y1+0.87*(y2-y1),(r'r$=%s$' % str(np.round(r,2))) +'\n'+('p$=%s$' % str(np.round(p,2))),color='red',va='center',ha='left')
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.savefig('size_vs_derrida_N%i.pdf' % len(Fs),bbox_inches = "tight")   










## FFLs

def compute_all_FFLs(Fs,Is,degrees,constantss,max_degree,N):
    import itertools
    all_max_degrees = list(map(max,degrees)) 
    
    ## find all FFLs
    all_ffls = []
    for i in range(len(Fs)):
        if all_max_degrees[i]>max_degree:
            all_ffls.append([])
            continue
        F = Fs[i]
        I = Is[i]
        A = can.adjacency_matrix(I,constantss[i])
        (ffls,types) = can.get_ffls(A,F,I)
        all_ffls.append(list(map(can.get_ffl_type_number,types)))
    
    all_ffls_flat = []
    for el in all_ffls:
        all_ffls_flat.extend(el)
        
    LEGEND = list(itertools.product(['decreasing', 'increasing'], repeat=3))
    LEGEND_COH = np.array(list(map(can.is_ffl_coherent,LEGEND)))
    
    nr_ffls = np.array([len(el) for el in all_ffls])
    nr_coherent_ffls = np.zeros(N,dtype=int)
    nr_incoherent_ffls = np.zeros(N,dtype=int)
    nr_unknown_ffls = np.zeros(N,dtype=int)
    nr_notreal_ffls = np.zeros(N,dtype=int)
    nr_specific_ffls = np.zeros((8,N),dtype=int)
    for ii,ffls in enumerate(all_ffls):
        for el in ffls:
            if el==-1:
                nr_unknown_ffls[ii]+=1
            elif el==-2:
                nr_notreal_ffls[ii]+=1
            elif LEGEND_COH[el] == True:
                nr_coherent_ffls[ii] += 1
                nr_specific_ffls[el,ii] += 1
            else:
                nr_incoherent_ffls[ii] += 1
                nr_specific_ffls[el,ii] += 1
    nr_real_ffls = nr_ffls-nr_notreal_ffls
    return nr_ffls,nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls,nr_specific_ffls,nr_real_ffls,nr_notreal_ffls

def plot_number_of_coherent_vs_incoherent_ffls_by_model(nr_real_ffls,nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls,N):
    #count coherent vs incoherent vs unknown stratified by model
    color_unknown = [0.5,0.5,0.5]
    color_incoh = [0.7,0.7,1]
    color_coh = [1,0.7,0.7]
    
    DONT_SHOW_ZERO_FFLS_NETWORKS = True
    LOG=False
    sorted_sizes = np.array(sorted(zip(nr_real_ffls,nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls),key=lambda x: (x[0],x[1]))).T
    if DONT_SHOW_ZERO_FFLS_NETWORKS:
        index = list(sorted_sizes[0]>0).index(True)
        sorted_sizes = sorted_sizes[:,index:]
    f,ax=plt.subplots(figsize=(10,4))
    ind = np.arange(sorted_sizes.shape[1])
    ax.bar(ind,sorted_sizes[1],width=0.7,color=color_coh,label='coherent FFLs',log=LOG)
    ax.bar(ind,sorted_sizes[2],width=0.7,bottom=sorted_sizes[1],color=color_incoh,label='incoherent FFLs',log=LOG)
    ax.bar(ind,sorted_sizes[3],width=0.7,bottom=sorted_sizes[1]+sorted_sizes[2],color=color_unknown,label='unknown type FFL',log=LOG)
    ax.set_xlim([-1,sorted_sizes.shape[1]])
    ax.set_ylabel('Number of FFLs')
    ax.set_xticks([])
    ax.set_xlabel('Gene regulatory networks')
    ax.legend(loc=2,frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig('number_of_coherent_ffls_N%i.pdf' % N)

def plot_proportion_of_coherent_vs_incoherent_ffls_by_model(nr_ffls,nr_real_ffls,nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls,N):
    #proportion coherent vs incoherent vs unknown stratified by model
    color_unknown = [0.5,0.5,0.5]
    color_incoh = [0.7,0.7,1]
    color_coh = [1,0.7,0.7]

    DONT_SHOW_ZERO_FFLS_NETWORKS = True
    LOG=False
    sorted_sizes = np.array(sorted(zip(nr_real_ffls,nr_coherent_ffls/nr_real_ffls,nr_incoherent_ffls/nr_real_ffls,nr_unknown_ffls/nr_real_ffls),key=lambda x: (x[0],x[1]))).T
    if DONT_SHOW_ZERO_FFLS_NETWORKS:
        index = list(sorted_sizes[0]>0).index(True)
        sorted_sizes = sorted_sizes[:,index:]
    f,ax=plt.subplots(figsize=(10,4))
    ind = np.arange(sorted_sizes.shape[1])
    ax.bar(ind,sorted_sizes[1],width=0.7,color=color_coh,label='coherent FFLs',log=LOG)
    ax.bar(ind,sorted_sizes[2],width=0.7,bottom=sorted_sizes[1],color=color_incoh,label='incoherent FFLs',log=LOG)
    ax.bar(ind,sorted_sizes[3],width=0.7,bottom=sorted_sizes[1]+sorted_sizes[2],color=color_unknown,label='unknown type FFL',log=LOG)
    ax.set_xlim([-1,sorted_sizes.shape[1]])
    ax.set_ylabel('Proportion of FFLs')
    ax.set_xticks([])
    ax.set_xlabel('Gene regulatory networks')
    ax2 = ax.twinx()
    ax2color = 'k'#[0.2,0.9,0.2]
    ax2.semilogy(ind,sorted_sizes[0],'-',lw=2,color=ax2color)
    ax2.set_ylabel('Number of FFLs',color=ax2color)
    ax2.tick_params(axis='y', labelcolor=ax2color)
    ax.legend(loc=8)
    plt.savefig('proportion_of_coherent_ffls_N%i.pdf' % N)

def plot_total_count_of_ffls(nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls,N):
    #total count, coherent vs incoherent vs unknown
    color_unknown = [0.5,0.5,0.5]
    color_incoh = [0.7,0.7,1]
    color_coh = [1,0.7,0.7]
    
    f,ax = plt.subplots()
    ax.bar([0],[sum(nr_coherent_ffls)],color=color_coh)
    ax.bar([1],[sum(nr_incoherent_ffls)],color=color_incoh)
    ax.bar([2],[sum(nr_unknown_ffls)],color=color_unknown)
    ax.set_xticks(range(3))
    ax.set_xticklabels(['coherent','incoherent','conditional'],rotation=90)
    ax.set_ylabel('total count')
    plt.savefig('total_count_of_ffls_N%i.pdf' % N)

def plot_total_count_of_specific_ffls_old(nr_specific_ffls,nr_unknown_ffls,N,figsize=(3.43,3.1)):
    import matplotlib
    import itertools
    def arrow_new(self, x, y, dx, dy, **kwargs):
        kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
        kwargs.setdefault('fc', 'black')
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)
        posA = x, y
        posB = x+dx, y+dy
        a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
        self.add_artist(a)
        return a

    #total count, specific FFL type
    LEGEND = list(itertools.product(['decreasing', 'increasing'], repeat=3))
    LEGEND_COH = np.array(list(map(can.is_ffl_coherent,LEGEND)))
    order = np.append(np.arange(8)[LEGEND_COH],np.arange(8)[~LEGEND_COH])
    order=np.array([7,4,2,1,5,6,3,0])
    cmap = matplotlib.cm.Paired
    colors = [cmap(i) for i in [3,2,1,0,5,4,7,6]]
    color_unknown = [0.5,0.5,0.5]
    
    f,ax = plt.subplots(figsize=figsize)
    width=0.8
    for i in range(8):
        ax.bar([i],[sum(nr_specific_ffls[order[i],:])],color=colors[i],width=width)
    ax.bar([8],[sum(nr_unknown_ffls)],color=color_unknown)
    ax.set_xticks(range(9))
    ax.set_xticklabels(['' for i in range(8)]+['conditional'],rotation=90)
    ax.xaxis.set_ticks_position('none') 
    ax.set_ylabel('Total count')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2 = ax.twinx()
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    yoffset = 0.12
    epsilon = 0.25
    total_length_y = 0.2
    activation_head_length = 5
    inhibition_head_length = 0.1
    ax2.plot([0-width/2,3+width/2],[1.15*(-yoffset-total_length_y),1.15*(-yoffset-total_length_y)],'k-',clip_on=False)
    ax2.plot([4-width/2,7+width/2],[1.15*(-yoffset-total_length_y),1.15*(-yoffset-total_length_y)],'k-',clip_on=False)
    ax2.text(1.5,1.35*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False,fontsize=9)
    ax2.text(5.6,1.35*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False,fontsize=9)
    ax2.set_ylim([0,1])
    for j in range(8):
        ax2.text(j,-.07,str(j+1),va='center',ha='center',color=colors[j],clip_on=False)
    for i in range(8):
        color=colors[i]
        direct,indirect1,indirect2 = LEGEND[order[i]]
        head_width = 4
        head_length = activation_head_length if direct == 'increasing' else inhibition_head_length
        arrow_new(ax2,i+epsilon,-yoffset,0,-total_length_y,clip_on=False,fc=color,ec=color,color=color,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        head_length = activation_head_length if indirect1 == 'increasing' else inhibition_head_length
        arrow_new(ax2,i+epsilon/4,-yoffset,-epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        head_length = activation_head_length if indirect2 == 'increasing' else inhibition_head_length
        arrow_new(ax2,i-epsilon,-yoffset-0.52*total_length_y,epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    ax2.set_yticks([])
    plt.gcf().subplots_adjust(bottom=0.3,left=0.2)
    plt.savefig('total_count_of_specific_ffls_nice_N%s.pdf' % str(N),bbox_inches = "tight")
    return order

def plot_proportion_of_specific_ffls_per_model(nr_real_ffls,nr_specific_ffls,nr_unknown_ffls,N):
    import matplotlib
    import itertools
    def arrow_new(self, x, y, dx, dy, **kwargs):
        kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
        kwargs.setdefault('fc', 'black')
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)
        posA = x, y
        posB = x+dx, y+dy
        a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
        self.add_artist(a)
        return a

    LEGEND = list(itertools.product(['decreasing', 'increasing'], repeat=3))
    LEGEND_COH = np.array(list(map(can.is_ffl_coherent,LEGEND)))
    order = np.append(np.arange(8)[LEGEND_COH],np.arange(8)[~LEGEND_COH])
    order=np.array([7,4,2,1,5,6,3,0])
    cmap = matplotlib.cm.Paired
    colors = [cmap(i) for i in [3,2,1,0,5,4,7,6]]
    color_unknown = [0.5,0.5,0.5]
    
    #proportion, specific FFL type stratified by model - sth odd is going on, not always adding up to 1 (unknown types?)
    DONT_SHOW_ZERO_FFLS_NETWORKS = True
    LOG=False
    cmap = matplotlib.cm.tab20c
    cmap = matplotlib.cm.Paired
    #cmap = matplotlib.cm.tab10
    sorted_sizes = np.array(sorted(zip(nr_real_ffls,nr_specific_ffls[0]/nr_real_ffls,nr_specific_ffls[1]/nr_real_ffls,nr_specific_ffls[2]/nr_real_ffls,nr_specific_ffls[3]/nr_real_ffls,nr_specific_ffls[4]/nr_real_ffls,nr_specific_ffls[5]/nr_real_ffls,nr_specific_ffls[6]/nr_real_ffls,nr_specific_ffls[7]/nr_real_ffls,nr_unknown_ffls/nr_real_ffls),key=lambda x: (x[0],x[1]))).T
    if DONT_SHOW_ZERO_FFLS_NETWORKS:
        index = list(sorted_sizes[0]>0).index(True)
        sorted_sizes = sorted_sizes[:,index:]

    #proportion, specific FFL type stratified by model - nice - sth odd is going on, not always adding up to 1 (unknown types?)
    DONT_SHOW_ZERO_FFLS_NETWORKS = True
    SHOW_LEGEND=False
    LOG=False
    cmap = matplotlib.cm.Paired
    color_ax2 = 'k'
    sorted_sizes = np.array(sorted(zip(nr_real_ffls,nr_specific_ffls[0]/nr_real_ffls,nr_specific_ffls[1]/nr_real_ffls,nr_specific_ffls[2]/nr_real_ffls,nr_specific_ffls[3]/nr_real_ffls,nr_specific_ffls[4]/nr_real_ffls,nr_specific_ffls[5]/nr_real_ffls,nr_specific_ffls[6]/nr_real_ffls,nr_specific_ffls[7]/nr_real_ffls,nr_unknown_ffls/nr_real_ffls),key=lambda x: (x[0],x[1]))).T
    if DONT_SHOW_ZERO_FFLS_NETWORKS:
        index = list(sorted_sizes[0]>0).index(True)
        sorted_sizes = sorted_sizes[:,index:]
    f,ax=plt.subplots(figsize=(7.3,4 if SHOW_LEGEND else 2.31))
    ind = np.arange(sorted_sizes.shape[1])
    bottom = np.zeros(sorted_sizes[order[0]].shape)
    for i in range(8):
        ax.bar(ind,sorted_sizes[1+order[i]],bottom=bottom,color=colors[i],log=LOG,width=0.7)
        bottom +=sorted_sizes[1+order[i]]     
    kk=9
    ax.bar(ind,sorted_sizes[kk],width=0.7,bottom=sum(sorted_sizes[1:kk]),color=color_unknown)
    ax.set_xlim([-1,sorted_sizes.shape[1]])
    ax.set_ylabel('Proportion of specific FFLs')
    ax.set_xticks([])
    ax.set_xlabel('Networks sorted by number of FFLs')
    ax2 = ax.twinx()
    ax2.semilogy(ind,sorted_sizes[0],'-',lw=3,color=color_ax2)
    ax2.set_ylabel('Number of FFLs',color=color_ax2)
    ax2.tick_params(axis='y', labelcolor=color_ax2)
    if SHOW_LEGEND:
        yoffset = -1.26
        epsilon = 0.25*((sorted_sizes.shape[1]+1)/9)
        total_length_y = 0.2
        activation_head_length = 7
        inhibition_head_length = 0.1
        ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
        ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
        ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
        ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
        ax.set_ylim([0,1])
        order = np.append(np.arange(8)[LEGEND_COH],np.arange(8)[~LEGEND_COH])
        for i,ii in enumerate(np.linspace(-1,sorted_sizes.shape[1],10)[1:-1]):
            direct,indirect1,indirect2 = LEGEND[order[i]]
            head_width = 7
            color=colors[i]
            head_length = activation_head_length if direct == 'increasing' else inhibition_head_length
            arrow_new(ax,ii+epsilon,-yoffset,0,-total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
            head_length = activation_head_length if indirect1 == 'increasing' else inhibition_head_length
            arrow_new(ax,ii+epsilon/4,-yoffset,-epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
            head_length = activation_head_length if indirect2 == 'increasing' else inhibition_head_length
            arrow_new(ax,ii-epsilon,-yoffset-0.52*total_length_y,epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        plt.gcf().subplots_adjust(top=0.75)
    plt.savefig('proportion_of_specific_type_ffls_nice_N%i.pdf' % N,bbox_inches = "tight")
    
    
    
def get_null_expectations_ffl(nr_specific_ffls,prop_pos_separate,N):
    import itertools
    LEGEND = list(itertools.product([0,1], repeat=3))
    nr_pos_per_type = np.sum(LEGEND,1)
    expected_number_per_type = np.zeros((N,8))
    for i in range(N):
        expected_number_per_type[i] = np.sum(nr_specific_ffls[:,i]) * prop_pos_separate[i]**nr_pos_per_type * (1-prop_pos_separate[i])**(3-nr_pos_per_type)
    return expected_number_per_type

def plot_total_count_of_specific_ffls(nr_specific_ffls,nr_unknown_ffls,expected_number_per_type,N,figsize=(3.43,3.1)):
    import matplotlib
    import itertools
    def arrow_new(self, x, y, dx, dy, **kwargs):
        kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
        kwargs.setdefault('fc', 'black')
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)
        posA = x, y
        posB = x+dx, y+dy
        a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
        self.add_artist(a)
        return a

    #total count, specific FFL type
    LEGEND = list(itertools.product(['decreasing', 'increasing'], repeat=3))
    LEGEND_COH = np.array(list(map(can.is_ffl_coherent,LEGEND)))
    order = np.append(np.arange(8)[LEGEND_COH],np.arange(8)[~LEGEND_COH])
    order=np.array([7,4,2,1,5,6,3,0])
    cmap = matplotlib.cm.Paired
    colors = [cmap(i) for i in [3,2,1,0,5,4,7,6]]
    color_unknown = [0.5,0.5,0.5]
    
    f,ax = plt.subplots(figsize=figsize)
    width=0.8
    for i in range(8):
        ax.bar([i],[sum(nr_specific_ffls[order[i],:])],color=colors[i],width=width)
        ax.plot([i-width/2+0.05,i+width/2-0.05],[sum(expected_number_per_type[:,order[i]]),sum(expected_number_per_type[:,order[i]])],'k-',label = 'expected' if i==0 else None)
    ax.bar([8],[sum(nr_unknown_ffls)],color=color_unknown)
    ax.legend(loc='best',frameon=False)
    ax.set_xticks(range(9))
    ax.set_xticklabels(['' for i in range(8)]+['conditional'],rotation=90)
    ax.xaxis.set_ticks_position('none') 
    ax.set_ylabel('Total count')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2 = ax.twinx()
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    yoffset = 0.12
    epsilon = 0.25
    total_length_y = 0.2
    activation_head_length = 5
    inhibition_head_length = 0.1
    ax2.plot([0-width/2,3+width/2],[1.15*(-yoffset-total_length_y),1.15*(-yoffset-total_length_y)],'k-',clip_on=False)
    ax2.plot([4-width/2,7+width/2],[1.15*(-yoffset-total_length_y),1.15*(-yoffset-total_length_y)],'k-',clip_on=False)
    ax2.text(1.5,1.35*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False,fontsize=9)
    ax2.text(5.6,1.35*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False,fontsize=9)
    ax2.set_ylim([0,1])
    for j in range(8):
        ax2.text(j,-.07,str(j+1),va='center',ha='center',color=colors[j],clip_on=False)
    for i in range(8):
        color=colors[i]
        direct,indirect1,indirect2 = LEGEND[order[i]]
        head_width = 4
        head_length = activation_head_length if direct == 'increasing' else inhibition_head_length
        arrow_new(ax2,i+epsilon,-yoffset,0,-total_length_y,clip_on=False,fc=color,ec=color,color=color,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        head_length = activation_head_length if indirect1 == 'increasing' else inhibition_head_length
        arrow_new(ax2,i+epsilon/4,-yoffset,-epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        head_length = activation_head_length if indirect2 == 'increasing' else inhibition_head_length
        arrow_new(ax2,i-epsilon,-yoffset-0.52*total_length_y,epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    ax2.set_yticks([])
    plt.gcf().subplots_adjust(bottom=0.3,left=0.2)
    plt.savefig('total_count_of_specific_ffls_nice_N%s.pdf' % str(N),bbox_inches = "tight")
    return order
















def compute_layer_regulator_vs_layer_intermediate_in_ffls(Fs,Is,constantss):
    ## find all FFLs and determine if the direct or the indirect regulator of the common target is more important in terms of canalization
    def build_dummy_vector(kis):
        res = []
        for i,el in enumerate(kis):
            res.extend([i+1]*el)
        return res
    
    res = []
    res_detailed = []
    nr_layers = 5 #set this or run once and then see what the max is
    res_mat = np.zeros((nr_layers+1,nr_layers+1),dtype=int)
    for i in range(len(Fs)):
        F = Fs[i]
        I = Is[i]
        A = can.adjacency_matrix(I,constantss[i])
        try:
            (ffls,types) = can.get_ffls(A,F,I)
        except OverflowError: #happens if F == [], i.e. if degree > max_degree
            continue #just ignore those few networks
        #res.append([])
        for [regulator,intermediate,target] in ffls:
            f = F[target]
            index_regulator = list(I[target]).index(regulator)
            index_intermediate = list(I[target]).index(intermediate)
            (n,k,can_inputs,can_outputs,corefunction,can_order) = can.get_canalizing_depth_inputs_outputs_corefunction_order(f)
            kis = can.get_layerstructure_given_canalizing_outputs_and_corefunction(can_outputs,corefunction,n)
            dummy_kis = build_dummy_vector(kis)
    
            try:
                index_in_canalizing_order_regulator = list(can_order).index(index_regulator)
                layer_regulator = dummy_kis[index_in_canalizing_order_regulator]
            except ValueError:
                layer_regulator = 1000
            try:
                index_in_canalizing_order_intermediate = list(can_order).index(index_intermediate)
                layer_intermediate = dummy_kis[index_in_canalizing_order_intermediate]
            except ValueError:
                layer_intermediate = 2000
            layers = [layer_regulator-1 if layer_regulator<1000 else nr_layers,layer_intermediate-1 if layer_intermediate<2000 else nr_layers]
            res_mat[layers[0],layers[1]] += 1
            res_detailed.append(','.join(list(map(str,layers))))
            res.append(layer_regulator-layer_intermediate)
    res = np.array(res)
    res_detailed = np.array(res_detailed)
    
    pd.DataFrame(res_mat,columns=['Layer %i' % (i+1) for i in range(nr_layers)]+['not canalizing'],index=['Layer %i' % (i+1) for i in range(nr_layers)]+['not canalizing']).to_excel('layer_regulator_row_vs_layer_intermediate_col_in_ffls_N%i.xlsx' % len(Fs))


def compute_effectiveness_regulator_vs_effectiveness_intermediate_in_ffls(Fs,Is,constantss,upto=8):
    import cana
    ## find all FFLs and determine if the direct or the indirect regulator of the common target is more important in terms of canalization
    
    res = []
    skipped_ffls = np.zeros(len(Fs))
    for i in range(len(Fs)):
        F = Fs[i]
        I = Is[i]
        A = can.adjacency_matrix(I,constantss[i])
        ffls = can.get_ffls(A)
        #res.append([])
        for [regulator,intermediate,target] in ffls:
            f = F[target]
            n = len(I[target])
            if n > upto: #happens if degree(f) > upto
                skipped_ffls[i]+=1
                continue
            
            index_regulator = list(I[target]).index(regulator)
            index_intermediate = list(I[target]).index(intermediate)
            edge_effectiveness = cana.boolean_node.BooleanNode(k=n,inputs=range(n),outputs=f).edge_effectiveness()
            edge_effectiveness_regulator = edge_effectiveness[index_regulator]
            edge_effectiveness_intermediate = edge_effectiveness[index_intermediate]
            
            effective_degree = sum([el>0 for el in edge_effectiveness])
            
            res.append([i,effective_degree,n,edge_effectiveness_regulator,edge_effectiveness_intermediate])
    res = np.array(res)
    return res,skipped_ffls
    
def plot_effectiveness_regulator_vs_effectiveness_intermediate_in_ffls(res,N,upto=8,figsize=(5,3.1)):
    import matplotlib.patches as mpatches
    import scipy.stats as stats
    labels = []
    def add_label(violin,label):
        color = violin['bodies'][0].get_facecolor().flatten()
        labels.append((mpatches.Patch(color = color),label))
        
    xs = np.arange(2,upto+1)
    f,ax = plt.subplots(figsize=figsize)
    add_label(ax.violinplot([res[res[:,1]==n,3] for n in xs],positions=xs-0.2,showextrema=False,showmeans=True,widths=0.35),'master regulator')
    add_label(ax.violinplot([res[res[:,1]==n,4] for n in xs],positions=xs+0.2,showextrema=False,showmeans=True,widths=0.35),'intermediate')
    for n in xs:
        wilcoxon_p = stats.wilcoxon(res[res[:,1]==n,3],res[res[:,1]==n,4])[1]
        print(n,wilcoxon_p)
        ax.text(n,1.07,'n = '+str(sum(res[:,1]==n)),va='center',ha='center',fontsize=7)
        ax.text(n,1,'p > 0.05' if wilcoxon_p>0.05 else 'p = '+"{:.0e}".format(wilcoxon_p),va='center',ha='center',fontsize=7)
    ax.set_xticks(xs)
    ax.set_ylim([-0.05,1.1])
    ax.legend(*zip(*labels),ncol=1,loc=3,frameon=False)#,bbox_to_anchor=[0.,-0.18])
    ax.set_xlabel("Essential in-degree of target gene")
    ax.set_ylabel("Edge effectiveness")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig('effectiveness_regulator_vs_effectiveness_intermediate_in_fflsp_N%i.pdf' % N,bbox_inches = "tight")










## Feedbackloops and triads

def compute_all_FBLs(Fs,Is,degrees,constantss,max_degree,max_loop=6):
    import networkx as nx
    all_max_degrees = list(map(max,degrees)) 

    all_loops = []
    all_types = []
    for i in range(len(Fs)):
        if all_max_degrees[i]>max_degree:
            all_types.append([])
            all_loops.append([])
            continue
        F = Fs[i]
        I = Is[i]
        edges = []
        for j,regulators in enumerate(I):
            if j>=n_variables[i]: #exclude constant self-loops
                break
            for ii in regulators:
                edges.append((ii,j))
        G=nx.DiGraph(edges)
        #triads.append(nx.triadic_census(G))
        loops = list(can.simple_cycles(G,max_loop))
        all_types.append([])
        for loop in loops:
            all_types[-1].append(can.get_type_of_loop(loop,F,I))
        all_loops.append(loops)
        
    #all_triads_keys = triads[0].keys()
    #all_triads_counts = [sum([el[key] for el in triads]) for key in all_triads_keys]
    
    all_loops_flat = []
    for el in all_loops:
        all_loops_flat.extend(el)

    nr_pos_loops = np.zeros((max_loop,N),dtype=int)
    nr_neg_loops = np.zeros((max_loop,N),dtype=int)
    nr_unknown_loops = np.zeros((max_loop,N),dtype=int)
    nr_notreal_loops = np.zeros((max_loop,N),dtype=int)
    nr_specific_k_loops = np.zeros((max_loop,max_loop+1,N),dtype=int)
    for ii,types in enumerate(all_types):
        for type_ in types:
            k = len(type_)
            el = can.get_loop_type_number(type_)
            if el==-1:
                nr_unknown_loops[k-1,ii]+=1
            elif el==-2:
                nr_notreal_loops[k-1,ii]+=1
            elif el%2==0:
                nr_pos_loops[k-1,ii] += 1
                nr_specific_k_loops[k-1,el,ii] += 1
            else:
                nr_neg_loops[k-1,ii] += 1
                nr_specific_k_loops[k-1,el,ii] += 1
            
    nr_real_loops = nr_pos_loops + nr_neg_loops + nr_unknown_loops
    nr_loops = nr_real_loops + nr_notreal_loops
    return nr_loops,nr_pos_loops,nr_neg_loops,nr_unknown_loops,nr_notreal_loops,nr_specific_k_loops,nr_real_loops,all_types,all_loops

def plot_number_of_pos_neg_and_unknown_FBLs(nr_real_loops,nr_pos_loops,nr_neg_loops,nr_unknown_loops,N,max_loop=6):
    #count positive vs negative vs unknown loops, stratified by model
    color_unknown = [0.5,0.5,0.5]
    color_neg = [0.7,0.7,1]
    color_pos = [1,0.7,0.7]
    
    for k in range(max_loop):
        DONT_SHOW_ZERO_LOOP_NETWORKS = True
        LOG=False
        sorted_sizes = np.array(sorted(zip(nr_real_loops[k],nr_pos_loops[k],nr_neg_loops[k],nr_unknown_loops[k]),key=lambda x: (x[0],x[1]))).T
        if DONT_SHOW_ZERO_LOOP_NETWORKS:
            index = list(sorted_sizes[0]>0).index(True)
            sorted_sizes = sorted_sizes[:,index:]
        f,ax=plt.subplots(figsize=(10,4))
        ind = np.arange(sorted_sizes.shape[1])
        ax.bar(ind,sorted_sizes[1],width=0.7,color=color_pos,label='positive %i-loops' % (k+1),log=LOG)
        ax.bar(ind,sorted_sizes[2],width=0.7,bottom=sorted_sizes[1],color=color_neg,label='negative %i-loops' % (k+1),log=LOG)
        ax.bar(ind,sorted_sizes[3],width=0.7,bottom=sorted_sizes[1]+sorted_sizes[2],color=color_unknown,label='unknown type %i-loops' % (k+1),log=LOG)
        ax.set_xlim([-1,sorted_sizes.shape[1]])
        ax.set_ylabel('Number of %i-loops' % (k+1))
        ax.set_xticks([])
        ax.set_xlabel('Gene regulatory networks')
        ax.legend(loc=2)
        plt.savefig('number_of_pos_neg_unknown_%iloops_N%i.pdf' % (k+1,N))
    
def plot_total_count_of_FBLs(nr_pos_loops,nr_neg_loops,nr_unknown_loops,N,max_loop=6):
    #total count positive vs negative vs unknown loops
    color_unknown = [0.5,0.5,0.5]
    color_neg = [0.7,0.7,1]
    color_pos = [1,0.7,0.7]
    for k in range(max_loop):
        f,ax = plt.subplots()
        ax.bar([0],[sum(nr_pos_loops[k])],color=color_pos)
        ax.bar([1],[sum(nr_neg_loops[k])],color=color_neg)
        ax.bar([2],[sum(nr_unknown_loops[k])],color=color_unknown)
        ax.set_xticks(range(3))
        ax.set_xticklabels(['positive','negative','conditional'],rotation=90)
        ax.set_ylabel('total number of %i-loops' % (k+1))
        plt.tight_layout()
        plt.savefig('total_count_of_%iloops_N%i.pdf' % (k+1,N))
    
def plot_total_count_of_specific_FBLs(nr_specific_k_loops,nr_unknown_loops,N,max_loop=6,LEGEND=False):
    import matplotlib
    color_unknown = [0.5,0.5,0.5]
    #total count, specific type of loops
    cmap = matplotlib.cm.tab20c
    for k in range(max_loop): 
        f,ax = plt.subplots(figsize=(1.3+0.3*k,2))
        order = np.array([0,2,4,6,1,3,5])
        order = order[order<k+2]
        for i in order:
            ax.bar([i],[sum(nr_specific_k_loops[k,i,:])],color=cmap(4*(i%2)+i//2),label=str(i))
        ax.bar([k+2],[sum(nr_unknown_loops[k])],color=color_unknown,label='n.a.')
        ax.set_xticks(range(k+3))
        if LEGEND:
            ax.legend(loc='best',frameon=False,title='number of\ninhibitory edges',ncol=2,bbox_to_anchor=[1.05,0.5])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticklabels([str(k+1-i)+r'$+$'+'\n'+str(i)+r'$-$' for i in range(k+2)]+['n.a.'])
        #ax.set_xticklabels([(r'$+$'*(k+1-i))+(r'$-$'*(i)) for i in range(k+2)]+['unknown'],rotation=90)
        #if k==0:
        #    ax.set_ylabel('Total count')#'number of %i-loops' % (k+1))
        plt.tight_layout()
        plt.savefig('total_count_of_specific_%iloops_N%i.pdf' % (k+1,N),bbox_inches = "tight")

def plot_total_count_of_pos_vs_neg_FBLs(nr_specific_k_loops,nr_unknown_loops,prop_pos_separate,prop_pos_specific_scc_of_loop2,N,max_loop=6):
    import matplotlib
    import scipy.special
    #total count, specific type of loops
    cmap = matplotlib.cm.tab20c

    PROPORTION=True
    height_rectangle = 0.1
    f,ax = plt.subplots(figsize=(4.2,2))
    colors = [cmap(el) for el in [0,4,1,5,2,6,3]]
    NULLMODEL = 1
    width = 0.62
    epsilon = 0.25
    lw=1
    for k in range(max_loop): 
        ax.text(3*k,1.05,'n=%i' % int(np.sum(nr_specific_k_loops[k,:,:])),va='center',ha='center',clip_on=False )

        order = np.array([0,2,4,6,5,3,1])
        #order = np.array([0,1,2,3,4,5,6])
        order = order[order<=k+1]
        s=0
        expected = np.array([sum(np.sum(nr_specific_k_loops[k,:,:],0) * scipy.special.binom(k+1,j) * prop_pos_separate**(k+1-j) * (1-prop_pos_separate)**(j))/(1 if not PROPORTION else np.sum(nr_specific_k_loops[k,:,:])) for j in range(k+2)])
        for el in order:
            ax.bar([3*k-width-epsilon],[expected[el]],color=colors[el],bottom=[s],width=width)
            s+=expected[el]
        
        s=0
        expected = np.array(prop_pos_specific_scc_of_loop2[k]).mean(0)
        for el in order:
            ax.bar([3*k],[expected[el]],color=colors[el],bottom=[s],width=width)
            s+=expected[el]
            
        observed =  np.sum(nr_specific_k_loops[k,:,:],1) if not PROPORTION else np.sum(nr_specific_k_loops[k,:,:],1)/np.sum(nr_specific_k_loops[k,:,:])       
        s=0
        for el in order:
            ax.bar([3*k+width+epsilon],[observed[el]],color=colors[el],bottom=[s],width=width)
            s+=observed[el]
    
    for k in range(max_loop):
        for y in [0,1]:
            ax.plot([3*k-3*width/2-epsilon,3*k-3*width/2-epsilon+width],[y,y],'k--',lw=lw)
            ax.plot([3*k-width/2,3*k-width/2+width],[y,y],'k:',lw=lw)
            ax.plot([3*k+width/2+epsilon,3*k+width/2+epsilon+width],[y,y],'k-',lw=lw)
        for dx in [0,width]:
            ax.plot([3*k-3*width/2-epsilon+dx,3*k-3*width/2-epsilon+dx],[0,1],'k--',lw=lw)
            ax.plot([3*k-width/2+dx,3*k-width/2+dx],[0,1],'k:',lw=lw)
            ax.plot([3*k+width/2+epsilon+dx,3*k+width/2+epsilon+dx],[0,1],'k-',lw=lw)
            
    # ax.text(17.2,0.87,'number\nnegative\nedges',ha='left',va='center')
    # for j in range(max_loop+1):
    #     ax.add_patch(matplotlib.patches.Rectangle([17.5 if j%2==0 else 19.5,0.56-0.17*(j//2)],0.5,height_rectangle,color=colors[j],clip_on=False))
    #     ax.text(0.8+(17.5 if j%2==0 else 19.5),0.56-0.17*(j//2)+height_rectangle/2.5,str(j),ha='left',va='center',clip_on=False)
    
    x1,x2 = ax.get_xlim()

    
    ax.text(0.7,-0.35,'number of inhibitory edges',ha='center',va='center')
    for j in range(max_loop+1):
        ax.add_patch(matplotlib.patches.Rectangle([-4+2.5*(j//2),-0.55 if j%2==0 else -0.71],1,height_rectangle,color=colors[j],clip_on=False))
        ax.text(-2.7+2.5*(j//2),(-0.55 if j%2==0 else -0.71)+height_rectangle/2.5,str(j),ha='left',va='center',clip_on=False)
    
    for j,(ls,label) in enumerate(zip(['--',':','-'],['expected (model 1)','expected (model 2)','observed'])):
        ax.text(9,-0.38-0.16*j+height_rectangle/2.5,label,ha='left',va='center',clip_on=False)
        ax.plot([7.7,8.7],[-0.38-0.16*j,-0.38-0.16*j],'k',ls=ls,lw=1,clip_on=False)
        ax.plot([7.7,8.7],np.array([-0.38-0.16*j,-0.38-0.16*j])+height_rectangle,'k',ls=ls,lw=1,clip_on=False)
        ax.plot([7.7,7.7],np.array([-0.38-0.16*j,-0.38-0.16*j+height_rectangle]),'k',ls=ls,lw=1,clip_on=False)
        ax.plot([8.7,8.7],np.array([-0.38-0.16*j,-0.38-0.16*j+height_rectangle]),'k',ls=ls,lw=1,clip_on=False)
    
    
    
    ax.set_ylim([0,1])
    ax.set_xlim([x1,x2])
    
    ax.set_xticks(np.arange(0,3*k+1,3))
    #ax.legend(loc='best',frameon=False,title='negative edges',ncol=2,bbox_to_anchor=[1.05,0.5])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticklabels([str(k+1) for k in range(max_loop)])
    ax.xaxis.set_ticks_position('none')
    ax.set_xlabel('Length of FBL')
    ax.set_ylabel('Proportion')
    
    #plt.tight_layout()
    plt.savefig('total_count_of_specific_%iloops_N%s_NULLMODEL%i.pdf' % (k+1,str(N),NULLMODEL),bbox_inches = "tight")





def get_null_expectations_fbl(Fs,Is,variabless,all_types,all_loops):
    import networkx as nx
    import scipy.special
    def associate_number_monotonic(*args):
        res = []
        for el in args:
            if el == 'decreasing':
                res.append(1)
            elif el== 'increasing':
                res.append(0)
            elif el== 'not essential':
                res.append(-2)        
            else:
                res.append(-1)        
        return res
    
    N = len(Fs)
    
    res = []
    prop_pos_separate = []
    prop_pos_separate_by_scc = []
    number_edges_separate_by_scc = []
    sccs = []
    total_edges_separate = []
    for I,var in zip(Is,variabless):
        n_var = len(var)
        total_edges_separate.append(sum(map(len,I[:n_var])))
    
    all_dict_sccs = []
    for i,(F,I,var) in enumerate(zip(Fs,Is,variabless)): 
        n_var = len(var)
        
        G = can.generate_networkx_graph_from_edges(I[:n_var],n_var)
        sccs.append([list(scc) for scc in nx.strongly_connected_components(G)])
        n_sccs = len(sccs[-1])
        dict_sccs = {}
        for i in range(n_sccs):
            for el in sccs[-1][i]:
                dict_sccs.update({el:i})
        all_dict_sccs.append(dict_sccs)
        dummy = [[] for _ in range(n_sccs)]
        dummy2 = []
        for target,f in enumerate(F[:n_var]):
            if len(f)==0: #happens if actual degree_f > max_degree
                continue
            (NONDEGENERATED,monotonic) = can.is_monotonic(f,True)
            nrs = associate_number_monotonic(*monotonic)
            n_essential=sum([el!=-2 for el in nrs])
            n = len(nrs)
            scc_target = dict_sccs[target]
            for el,regulator in zip(nrs,I[target]):
                try:
                    scc_regulator = dict_sccs[regulator]
                except KeyError: #happens if regulator is an external parameter, thus clearly not in the same SCC as the target
                    continue
                if scc_target==scc_regulator and el in [0,1]:
                    dummy[scc_target].append(el)
            for el in nrs:
                res.append([i,n_essential,n,el])
                if el in [0,1]:
                    dummy2.append(el)            
        prop_pos_separate.append((len(dummy2)-sum(dummy2))/len(dummy2))
        prop_pos_separate_by_scc.append([(len(vec)-sum(vec))/len(vec) if vec!=[] else np.nan for vec in dummy])
        number_edges_separate_by_scc.append([len(vec) for vec in dummy])
    res = np.array(res)    
    prop_pos_separate = np.array(prop_pos_separate)
    prop_pos_separate_by_scc = np.array(prop_pos_separate_by_scc)
    number_edges_separate_by_scc = np.array(number_edges_separate_by_scc)
    total_edges_separate = np.array(total_edges_separate)
    prop_edges_in_largest_scc = np.array(list(map(max,number_edges_separate_by_scc)))/total_edges_separate
    prop_pos_in_largest_scc = [prop_pos_separate_by_scc[i][np.argmax(number_edges_separate_by_scc[i])] for i in range(N)]
    size_largest_scc = [len(sccs[i][np.argmax(number_edges_separate_by_scc[i])]) for i in range(N)]
    
    max_loop=6
    nr_specific_k_loops = np.zeros((max_loop,max_loop+1,N),dtype=int)
    prop_pos_specific_scc_of_loop = [[[[] for ii in range(N)] for j in range(max_loop+1)] for i in range(max_loop)]
    prop_pos_specific_scc_of_loop2 = [[] for i in range(max_loop)]
    for ii,(types,loops) in enumerate(zip(all_types,all_loops)):
        for (type_,loop) in zip(types,loops):
            k = len(type_)
            el = can.get_loop_type_number(type_)
            if el<0:
                continue
            nr_neg = el
            nr_pos = k-nr_neg
            nr_specific_k_loops[k-1,el,ii] += 1
            id_scc = all_dict_sccs[ii][loop[0]]
            p_pos = prop_pos_separate_by_scc[ii][id_scc]
            prop_pos_specific_scc_of_loop2[k-1].append( [scipy.special.binom(k,j) * p_pos**(k-j) * (1-p_pos)**j for j in range(k+1)] )
            
            prop_pos = p_pos**nr_pos * (1-p_pos)**nr_neg
            prop_neg = p_pos**nr_neg * (1-p_pos)**nr_pos
            prop_pos_specific_scc_of_loop[k-1][el][ii].append(prop_pos/(prop_pos+prop_neg))#p_pos)
        for i in range(max_loop):
            for j in range(max_loop+1):
                if len(prop_pos_specific_scc_of_loop[i][j][ii]) > 0:
                    prop_pos_specific_scc_of_loop[i][j][ii] = np.mean(prop_pos_specific_scc_of_loop[i][j][ii])
                else:
                    prop_pos_specific_scc_of_loop[i][j][ii] = np.nan
    prop_pos_specific_scc_of_loop = np.array(prop_pos_specific_scc_of_loop)
        
    A = pd.DataFrame(np.c_[prop_pos_separate,prop_pos_in_largest_scc,size_largest_scc,list(map(len,variabless)),prop_edges_in_largest_scc,np.sum(nr_specific_k_loops,1)[3],np.sum(nr_specific_k_loops,1)[5]])
    A.to_excel('out.xlsx')
    return prop_pos_separate,prop_pos_separate_by_scc,prop_pos_specific_scc_of_loop,prop_pos_specific_scc_of_loop2

def plot_pos_vs_neg_regulations_in_FBLs_best(Fs,nr_specific_k_loops,all_loops,prop_pos_separate,prop_pos_specific_scc_of_loop,N):
    import matplotlib
    #expected proportion based on proportion pos vs neg regulations

    
    def arrow_new(self, x, y, dx, dy, **kwargs):
        kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
        kwargs.setdefault('fc', 'black')
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)
        posA = x, y
        posB = x+dx, y+dy
        a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
        self.add_artist(a)
        return a

    
    prop_pos = sum(res[:,3]==0)/(sum(res[:,3]==0)+sum(res[:,3]==1))
    neg_of_ns = [[4,4],[3,4],[6,6],[5,6],[4,6]]
    
        
    f,ax = plt.subplots(figsize=(4.2,2))
    width=0.65
    epsilon = 0.22
    colors = ['blue','orange']
    
    cmap = matplotlib.cm.tab20c
    colors = [cmap(el) for el in [0,4,1,5,2,6,3]]
    
    for i,neg_of_n in enumerate(neg_of_ns):
        nr_neg,k = neg_of_n
        nr_pos = k-nr_neg
        more_pos = nr_specific_k_loops[k-1,nr_pos,:]
        more_neg = nr_specific_k_loops[k-1,nr_neg,:]
        total_nr_loops = more_pos+more_neg
        
        sum_more_pos = sum(more_pos)
        sum_more_neg = sum(more_neg)
        sum_total_nr_loops = sum(total_nr_loops)
        
        ax.bar([3*i+width+epsilon],[sum_more_pos/sum_total_nr_loops],width=width,color=colors[nr_pos])
        ax.bar([3*i+width+epsilon],[sum_more_neg/sum_total_nr_loops],bottom=[sum_more_pos/sum_total_nr_loops],width=width,color=colors[nr_neg])
        ax.text(3*i,1.05,'n=%i' % sum_total_nr_loops,va='center',ha='center',clip_on=False )
        #ax.text(3*i-width/2-epsilon,sum_more_pos/sum_total_nr_loops/2,str(nr_neg)+r'$+$'+'\n'+str(nr_pos)+r'$-$',va='center',ha='center')
        #ax.text(3*i-width/2-epsilon,sum_more_pos/sum_total_nr_loops + sum_more_neg/sum_total_nr_loops/2,str(nr_pos)+r'$+$'+'\n'+str(nr_neg)+r'$-$',va='center',ha='center')
        for iii,NULLMODEL in enumerate([2,3]):
            if NULLMODEL==1:
                expected_more_pos = prop_pos**nr_neg*(1-prop_pos)**nr_pos
                expected_more_neg = prop_pos**nr_pos*(1-prop_pos)**nr_neg
            elif NULLMODEL==2:
                expected_more_pos = np.dot(total_nr_loops,prop_pos_separate**nr_neg * (1-prop_pos_separate)**nr_pos)
                expected_more_neg = np.dot(total_nr_loops,prop_pos_separate**nr_pos * (1-prop_pos_separate)**nr_neg)
            elif NULLMODEL==3:
                dummy = prop_pos_specific_scc_of_loop[k-1,nr_pos]
                isnotnan = np.isnan(prop_pos_specific_scc_of_loop[k-1,nr_neg])==False
                dummy[isnotnan] = 1-prop_pos_specific_scc_of_loop[k-1,nr_neg][isnotnan]
                expected_more_pos = np.nansum(total_nr_loops * dummy)#prop_pos_specific_scc_of_loop[k-1,nr_pos]**nr_neg * (1-prop_pos_specific_scc_of_loop[k-1,nr_pos])**nr_pos)
                expected_more_neg = np.nansum(total_nr_loops * (1-dummy))#prop_pos_specific_scc_of_loop[k-1,nr_pos]**nr_pos * (1-prop_pos_specific_scc_of_loop[k-1,nr_pos])**nr_neg)
            total_expected = expected_more_pos+expected_more_neg
            #print(NULLMODEL,neg_of_n,expected_more_pos/total_expected*100,expected_more_neg/total_expected*100)
            ax.bar([3*i-width-epsilon+iii*(width+epsilon)],[expected_more_pos/total_expected],width=width,color=colors[nr_pos])
            ax.bar([3*i-width-epsilon+iii*(width+epsilon)],[expected_more_neg/total_expected],bottom=[expected_more_pos/total_expected],width=width,color=colors[nr_neg])
        

    lw=1
    for k,_ in enumerate(neg_of_ns):
        for y in [0,1]:
            ax.plot([3*k-3*width/2-epsilon,3*k-3*width/2-epsilon+width],[y,y],'k--',lw=lw)
            ax.plot([3*k-width/2,3*k-width/2+width],[y,y],'k:',lw=lw)
            ax.plot([3*k+width/2+epsilon,3*k+width/2+epsilon+width],[y,y],'k-',lw=lw)
        for dx in [0,width]:
            ax.plot([3*k-3*width/2-epsilon+dx,3*k-3*width/2-epsilon+dx],[0,1],'k--',lw=lw)
            ax.plot([3*k-width/2+dx,3*k-width/2+dx],[0,1],'k:',lw=lw)
            ax.plot([3*k+width/2+epsilon+dx,3*k+width/2+epsilon+dx],[0,1],'k-',lw=lw)
            
    ax.plot([-3*width/2-epsilon,3*1+3*width/2+epsilon],[-0.16,-0.16],'k-',clip_on=False,lw=0.5)
    ax.text(1.5,-0.23,'4-loops',ha='center',va='center')
    ax.plot([3*2-3*width/2-epsilon,3*4+3*width/2+epsilon],[-0.16,-0.16],'k-',clip_on=False,lw=0.5)
    ax.text(9,-0.23,'6-loops',ha='center',va='center')
    
    ax.set_ylim([0,1])
    #ax.set_xlim([-1,3*len(neg_of_ns)-1])
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('none')
    ax.xaxis.set_ticks(3*np.arange(len(neg_of_ns)))
    ax.xaxis.set_ticklabels(['positive' if el[0]%2==0 else 'negative' for el in neg_of_ns])
    ax.xaxis.set_ticklabels(['%i vs %i' % (el[0],el[1]-el[0]) for el in neg_of_ns])
    ax.set_ylabel('Proportion')
    plt.savefig('pos_vs_neg_regulations_in_loop_N%s_nice_NULLMODELall.pdf' % (str(N)),bbox_inches = "tight")



## FFL clusters
def compute_all_FFL_clusters(Is,constantss,n_variables,degrees,max_degree,N,DEBUG = False):
    all_max_degrees = list(map(max,degrees)) 

    type_ffl_ffls = []
    
    for ii in range(N):
        if all_max_degrees[ii]>max_degree:
            type_ffl_ffls.append([])
            continue
        # F = Fs[ii]
        # if len(F)>400:
        #     continue
        # I = Is[ii]
        A = can.adjacency_matrix(Is[ii],constantss[ii])
        #(ffls,types) = can.get_ffls(A,F,I)
        ffls = can.get_ffls(A)
        
        type_ffl_ffl = []
        #Determines types of motifs from 2018 Science Advances paper and new motifs (FFL and 2-loop)
        n_ffls = len(ffls)
        for i in range(n_ffls):
            for j in range(i+1,n_ffls):
                intersect = set.intersection(set(ffls[i]),set(ffls[j]))
                if len( intersect )>0:
                    if len(intersect)==3: #new types not included in 2018 paper
                        types = dict(zip(ffls[i],['r','i','t']))
                        for new_type,el in zip(['r','i','t'],ffls[j]):
                            try:
                                type_el = types[el]
                                if new_type != type_el:
                                    types.update({el:'m'})
                            except:
                                types.update({el:new_type})
                        sorted_types = list(types.values())
                        sorted_types.sort()
                        string_type = ''.join(sorted_types)
                        if DEBUG:
                            print(string_type)
                        if string_type == 'mmr':
                            type_ffl_ffl.append(13)
                        elif string_type == 'mmt':
                            type_ffl_ffl.append(14)
                        elif string_type == 'mmm': #When all the nodes have switched roles
                            if ffls[i][1] == ffls[j][0] or ffls[i][1] == ffls[j][2]:
                                type_ffl_ffl.append(15)
                            else:
                                type_ffl_ffl.append(17)
                    elif len(intersect)==2: #types 7-12 from 2018 paper
                        types = dict(zip(ffls[i],['r','i','t']))
                        for new_type,el in zip(['r','i','t'],ffls[j]):
                            try:
                                type_el = types[el]
                                if new_type != type_el:
                                    types.update({el:'m'})
                            except:
                                types.update({el:new_type})
                        sorted_types = list(types.values())
                        sorted_types.sort()
                        string_type = ''.join(sorted_types)
                        if DEBUG:
                            print(string_type)
                        if string_type == 'irrt':
                            type_ffl_ffl.append(7)
                        elif string_type == 'imrt':
                            if ffls[i][2] == ffls[j][2]:#if the target in both FFLs is the same
                                type_ffl_ffl.append(8)
                            else:
                                type_ffl_ffl.append(11)
                        elif string_type == 'mmrt':
                            type_ffl_ffl.append(9)
                        elif string_type == 'iirt':
                            type_ffl_ffl.append(10)                
                        elif string_type == 'irtt':
                            type_ffl_ffl.append(12)
                    elif len(intersect)==1: #types 1-6 from 2018 paper
                        types = dict(zip(ffls[i],['r','i','t']))
                        for new_type,el in zip(['r','i','t'],ffls[j]):
                            try:
                                type_el = types[el]
                                if new_type != type_el:
                                    types.update({el:'m'})
                            except:
                                types.update({el:new_type})
                        sorted_types = list(types.values())
                        sorted_types.sort()
                        string_type = ''.join(sorted_types)
                        if DEBUG:
                            print(string_type)
                        if string_type == 'iirrt':
                            type_ffl_ffl.append(1)
                        elif string_type == 'imrrt':
                            type_ffl_ffl.append(2)
                        elif string_type == 'iimrt':
                            type_ffl_ffl.append(3)
                        elif string_type == 'irrtt':
                            type_ffl_ffl.append(4)                
                        elif string_type == 'imrtt':
                            type_ffl_ffl.append(5)
                        elif string_type == 'iirtt':
                            type_ffl_ffl.append(6)                        
        type_ffl_ffls.append(type_ffl_ffl)
        
        # #FFLs clustered with 3 loops
        # max_loop_length = 3
        # G=can.generate_networkx_graph_from_edges(I,n_variables[ii])
        # loops = list(can.simple_cycles(G,max_len = max_loop_length))
        # two_loops = [el for el in loops if len(el)==2]
        
        # for i in range(n_ffls):
        #     for j in range(len(two_loops)):            
        #         intersect = set.intersection(set(ffls[i]),set(two_loops[j]))
        #         if len( intersect )>0:
        #             types = dict(zip(ffls[i],['r','i','t']))
        #             for el in list(intersect):
        #                 types.update({el:'m'})
        #             sorted_types = list(types.values())
        #             sorted_types.sort()
        #             string_type = ''.join(sorted_types) 
        #             #print(string_type)
        #             if string_type == 'mmr':
        #                 type_ffl_2l.append(13)
        #             elif string_type == 'mmt':
        #                 type_ffl_2l.append(14)
        #             elif string_type == 'imm':
        #                 type_ffl_2l.append(18)
        #             elif string_type == 'imt':
        #                 type_ffl_2l.append(19)
        #             elif string_type == 'mrt':
        #                 type_ffl_2l.append(20)
        #             elif string_type == 'imr':
        #                 type_ffl_2l.append(21)
    
    all_type_ffl_ffls = []
    for i,el in enumerate(type_ffl_ffls): 
        all_type_ffl_ffls.extend(el)
        
    all_types = list(set(all_type_ffl_ffls))
    n_types = len(all_types)
    dict_types = dict(zip(all_types,range(n_types)))
    
    nr_specific_ffl_ffl = np.zeros((n_types,N),dtype=int)
    for i,el in enumerate(type_ffl_ffls): 
        for typ in el:
            nr_specific_ffl_ffl[dict_types[typ],i]+=1
    nr_total_ffl_ffl = np.sum(nr_specific_ffl_ffl,0)
    return nr_total_ffl_ffl,nr_specific_ffl_ffl

def plot_total_count_of_specific_type_ffls_clusters(nr_total_ffl_ffl,nr_specific_ffl_ffl,N):
    #total count, specific FFL type - nice
    import matplotlib
    def arrow_new(self, x, y, dx, dy, **kwargs):
        kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
        kwargs.setdefault('fc', 'black')
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)
        posA = x, y
        posB = x+dx, y+dy
        a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
        self.add_artist(a)
        return a
    
    types_graphic = [list([[0,1],[0,2],[1,2],[3,2],[3,4],[4,2]]),
                     list([[0,1],[0,2],[2,1],[3,2],[3,4],[4,2]]),
                     list([[0,1],[0,2],[1,2],[2,3],[3,4],[2,4]]),
                     list([[0,1],[0,2],[2,1],[3,2],[3,4],[2,4]]),
                     list([[0,1],[2,0],[2,1],[3,2],[3,4],[2,4]]),
                     list([[0,1],[2,0],[2,1],[2,3],[3,4],[2,4]]),
                     list([[0,1],[0,2],[1,2],[3,1],[3,2]]),
                     list([[0,1],[0,2],[1,2],[1,3],[3,2]]),
                     list([[0,1],[0,2],[1,2],[1,3],[2,3]]),
                     list([[1,0],[0,2],[1,2],[1,3],[3,2]]),
                     list([[1,0],[2,0],[1,2],[1,3],[3,2]]),
                     list([[1,0],[2,0],[1,2],[1,3],[2,3]]),
                     list([[0,1],[0,2],[1,2],[2,1]]),
                     list([[0,1],[0,2],[1,2],[1,0]]),
                     list([[0,1],[0,2],[1,2],[1,0],[2,0]])]
    type_graphic = types_graphic[0]
    
    n_types=15
    order=list(range(n_types))
    
    cmap = matplotlib.cm.tab20
    
    COLOR_BY_TYPE = True
    colors_r_i_t_m = [cmap(0),cmap(2),cmap(5),[0.5,0.5,0.5]]
    
    #nice for publication
    f,ax = plt.subplots(figsize=(7,3.5))
    width=0.8
    for i in range(n_types):
        ax.bar([i],[sum(nr_specific_ffl_ffl[order[i],:])],color=cmap(i) if COLOR_BY_TYPE else 'k',width=width)
    ax.set_xticks(range(n_types))
    ax.set_xticklabels(['' for i in range(n_types)],rotation=90)
    ax.xaxis.set_ticks_position('none') 
    ax.set_ylabel('Total count')
    ax.set_xlim([-.5,n_types-.5])
    
    ax2 = ax.twinx()
    
    ax2.set_ylim([0,1])    
    delta_y = 0.13
    delta_x = 0.4*width
    #ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    #ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    #ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
    #ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
    #ax.set_ylim([0,1])
    
    head_width = 4
    head_length = 4
    
    ycenter = -0.31
    points = [(-1,1),(1,1),(0,0),(-1,-1),(1,-1)] #from top to bottom, left to right
    for j in range(n_types):
        ax2.text(j,-0.1,str(j+1),va='center',ha='center',color=cmap(j) if COLOR_BY_TYPE else 'k',clip_on=False)
    for j in range(6):
        xcenter=j
        color=cmap(j) if COLOR_BY_TYPE else 'k'
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            #ax.text(xcenter-0.8*delta_x,ycenter,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        for ii,point in enumerate(points):
            if not COLOR_BY_TYPE:
                ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
                if ii not in ffls_in_cluster[:,1:]: #regulator
                    color=colors_r_i_t_m[0]
                elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                    color=colors_r_i_t_m[1]
                elif ii not in ffls_in_cluster[:,:2]: #target
                    color=colors_r_i_t_m[2]
                else:
                    color=colors_r_i_t_m[3] 
            ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    points = [(0,1),(-1,0),(1,0),(0,-1)] #from top to bottom, left to right
    for j in range(6,12):
        xcenter=j
        color=cmap(j) if COLOR_BY_TYPE else 'k'
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        for ii,point in enumerate(points):
            if not COLOR_BY_TYPE:
                ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
                if ii not in ffls_in_cluster[:,1:]: #regulator
                    color=colors_r_i_t_m[0]
                elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                    color=colors_r_i_t_m[1]
                elif ii not in ffls_in_cluster[:,:2]: #target
                    color=colors_r_i_t_m[2]
                else:
                    color=colors_r_i_t_m[3]                
            ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    points = [(0,.5),(-1,-.5),(1,-.5)] #from top to bottom, left to right
    for j in range(12,15):
        xcenter=j
        color=cmap(j) if COLOR_BY_TYPE else 'k'
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        for ii,point in enumerate(points):
            if not COLOR_BY_TYPE:
                ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
                if ii not in ffls_in_cluster[:,1:]: #regulator
                    color=colors_r_i_t_m[0]
                elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                    color=colors_r_i_t_m[1]
                elif ii not in ffls_in_cluster[:,:2]: #target
                    color=colors_r_i_t_m[2]
                else:
                    color=colors_r_i_t_m[3] 
            ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    if not COLOR_BY_TYPE:
        positions = [0,3.6,6.8,9.8]
        for j,text in enumerate(['master regulator','intermediate','target','mixed']):
            ax2.text(positions[j]+j*0.5+0.3,-0.57,text,va='center',ha='left',color=cmap(j) if COLOR_BY_TYPE else 'k',clip_on=False)
            ax2.plot([positions[j]+j*0.5],[-0.57],'o',color=colors_r_i_t_m[j],markersize=7,clip_on=False)
    ax2.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    plt.gcf().subplots_adjust(bottom=0.3,left=0.12,right=0.99)
    plt.savefig('total_count_of_specific_ffls_ffls_nice_N%i.pdf' % (N),bbox_inches = "tight")
    
    
def plot_total_count_of_specific_type_ffls_clusters_firsttwelve(nr_total_ffl_ffl,nr_specific_ffl_ffl,N):
    #total count, specific FFL type - nice
    import matplotlib
    def arrow_new(self, x, y, dx, dy, **kwargs):
        kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
        kwargs.setdefault('fc', 'black')
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)
        posA = x, y
        posB = x+dx, y+dy
        a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
        self.add_artist(a)
        return a
    
    types_graphic = [list([[0,1],[0,2],[1,2],[3,2],[3,4],[4,2]]),
                     list([[0,1],[0,2],[2,1],[3,2],[3,4],[4,2]]),
                     list([[0,1],[0,2],[1,2],[2,3],[3,4],[2,4]]),
                     list([[0,1],[0,2],[2,1],[3,2],[3,4],[2,4]]),
                     list([[0,1],[2,0],[2,1],[3,2],[3,4],[2,4]]),
                     list([[0,1],[2,0],[2,1],[2,3],[3,4],[2,4]]),
                     list([[0,1],[0,2],[1,2],[3,1],[3,2]]),
                     list([[0,1],[0,2],[1,2],[1,3],[3,2]]),
                     list([[0,1],[0,2],[1,2],[1,3],[2,3]]),
                     list([[1,0],[0,2],[1,2],[1,3],[3,2]]),
                     list([[1,0],[2,0],[1,2],[1,3],[3,2]]),
                     list([[1,0],[2,0],[1,2],[1,3],[2,3]]),
                     list([[0,1],[0,2],[1,2],[2,1]]),
                     list([[0,1],[0,2],[1,2],[1,0]]),
                     list([[0,1],[0,2],[1,2],[1,0],[2,0]])]
    type_graphic = types_graphic[0]
    
    n_types=12
    order=list(range(n_types))
    
    cmap = matplotlib.cm.tab20
    
    COLOR_BY_TYPE = False
    colors_r_i_t_m = [cmap(0),cmap(2),cmap(5),[0.5,0.5,0.5]]
    
    #nice for publication
    f,ax = plt.subplots(figsize=(7,3.5))
    width=0.8
    for i in range(n_types):
        ax.bar([i],[sum(nr_specific_ffl_ffl[order[i],:])],color=cmap(i) if COLOR_BY_TYPE else 'k',width=width)
    ax.set_xticks(range(n_types))
    ax.set_xticklabels(['' for i in range(n_types)],rotation=90)
    ax.xaxis.set_ticks_position('none') 
    ax.set_ylabel('Total count')
    ax.set_xlim([-.5,n_types-.5])
    
    ax2 = ax.twinx()
    
    ax2.set_ylim([0,1])    
    delta_y = 0.13
    delta_x = 0.4*width
    #ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    #ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    #ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
    #ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
    #ax.set_ylim([0,1])
    
    head_width = 4
    head_length = 4
    
    ycenter = -0.31
    points = [(-1,1),(1,1),(0,0),(-1,-1),(1,-1)] #from top to bottom, left to right
    for j in range(12):
        ax2.text(j,-0.1,str(j+1),va='center',ha='center',color=cmap(j) if COLOR_BY_TYPE else 'k',clip_on=False)
    for j in range(6):
        xcenter=j
        color=cmap(j) if COLOR_BY_TYPE else 'k'
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            #ax.text(xcenter-0.8*delta_x,ycenter,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        for ii,point in enumerate(points):
            if not COLOR_BY_TYPE:
                ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
                if ii not in ffls_in_cluster[:,1:]: #regulator
                    color=colors_r_i_t_m[0]
                elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                    color=colors_r_i_t_m[1]
                elif ii not in ffls_in_cluster[:,:2]: #target
                    color=colors_r_i_t_m[2]
                else:
                    color=colors_r_i_t_m[3] 
            ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    points = [(0,1),(-1,0),(1,0),(0,-1)] #from top to bottom, left to right
    for j in range(6,12):
        xcenter=j
        color=cmap(j) if COLOR_BY_TYPE else 'k'
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        for ii,point in enumerate(points):
            if not COLOR_BY_TYPE:
                ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
                if ii not in ffls_in_cluster[:,1:]: #regulator
                    color=colors_r_i_t_m[0]
                elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                    color=colors_r_i_t_m[1]
                elif ii not in ffls_in_cluster[:,:2]: #target
                    color=colors_r_i_t_m[2]
                else:
                    color=colors_r_i_t_m[3]                
            ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    points = [(0,.5),(-1,-.5),(1,-.5)] #from top to bottom, left to right
    for j in range(12,15):
        xcenter=j
        color=cmap(j) if COLOR_BY_TYPE else 'k'
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        for ii,point in enumerate(points):
            if not COLOR_BY_TYPE:
                ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
                if ii not in ffls_in_cluster[:,1:]: #regulator
                    color=colors_r_i_t_m[0]
                elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                    color=colors_r_i_t_m[1]
                elif ii not in ffls_in_cluster[:,:2]: #target
                    color=colors_r_i_t_m[2]
                else:
                    color=colors_r_i_t_m[3] 
            ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    positions = [0,3.6,6.8,9.8]
    for j,text in enumerate(['master regulator','intermediate','target','mixed']):
        ax2.text(positions[j]+0.3,-0.57,text,va='center',ha='left',color=cmap(j) if COLOR_BY_TYPE else 'k',clip_on=False)
        ax2.plot([positions[j]],[-0.57],'o',color=colors_r_i_t_m[j],markersize=7,clip_on=False)
    ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    ax2.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    plt.gcf().subplots_adjust(bottom=0.3,left=0.12,right=0.99)
    plt.savefig('total_count_of_specific_ffls_ffls_nice_N%i_firsttwelve.pdf' % (N),bbox_inches = "tight")
    
    
    
    
def plot_total_count_of_specific_type_ffls_clusters_lastthree(nr_total_ffl_ffl,nr_specific_ffl_ffl,N):
    #total count, specific FFL type - nice
    import matplotlib
    
    n_types=3
    
    cmap = matplotlib.cm.tab20
    
    COLOR_BY_TYPE = False
    
    #nice for publication
    f,ax = plt.subplots(figsize=(1.7,1.35))
    width=0.8
    for i in range(n_types):
        ax.bar([i],[sum(nr_specific_ffl_ffl[12+i,:])],color=cmap(i) if COLOR_BY_TYPE else 'k',width=width)
    ax.set_xticks(range(n_types))
    ax.set_xticklabels(['' for i in range(n_types)],rotation=90)
    ax.xaxis.set_ticks_position('none') 
    ax2 = ax.twinx()
    ax2.set_ylim([0,1]) 
    for j in range(3):
        ax2.text(j,-0.2,str(12+j+1),va='center',ha='center',color=cmap(j) if COLOR_BY_TYPE else 'k',clip_on=False)
    ax.set_ylabel('Total count')
    ax.set_xlim([-.5,n_types-.5])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)    
    ax2.set_yticks([])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.gcf().subplots_adjust(bottom=0.3,left=0.12,right=0.99)
    plt.savefig('total_count_of_specific_ffls_ffls_nice_N%i_lastthree.pdf' % (N),bbox_inches = "tight")
    
    
    
    
    
    
def draw_ffl():
    #total count, specific FFL type - nice
    import matplotlib
    def arrow_new(self, x, y, dx, dy, **kwargs):
        kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
        kwargs.setdefault('fc', 'black')
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)
        posA = x, y
        posB = x+dx, y+dy
        a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
        self.add_artist(a)
        return a
    
    cmap = matplotlib.cm.tab20
    
    COLOR_BY_TYPE = False
    colors_r_i_t_m = [cmap(0),cmap(2),cmap(5),[0.5,0.5,0.5]]
    
    #nice for publication
    f,ax = plt.subplots(figsize=(7,3.5))
    ax2 = ax.twinx()
    width=0.8
    ax2.set_ylim([0,1])    
    delta_y = 0.13
    delta_x = 0.4*width
    #ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    #ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    #ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
    #ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
    #ax.set_ylim([0,1])
    
    head_width = 4
    head_length = 4
    
    ycenter = -0.31
    points = [(0,1),(-1.5,0),(0,-1)] #from top to bottom, left to right
    j=0
    xcenter=j
    color=cmap(j) if COLOR_BY_TYPE else 'k'
    
    types_graphic = [list([[0,1],[0,2],[1,2]])]*12
    for j in range(12):
        xcenter=j
        color=cmap(j) if COLOR_BY_TYPE else 'k'
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            color = colors_r_i_t_m[type_graphic[i][0]] if type_graphic[i][1]==2 else 'k'
            #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        for ii,point in enumerate(points):
            color=colors_r_i_t_m[ii] 
            ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    positions = [0,3.6,6.8,9.8]
    for j,text in enumerate(['master regulator','intermediate','target','mixed']):
        ax2.text(positions[j]+0.3,-0.57,text,va='center',ha='left',color=cmap(j) if COLOR_BY_TYPE else 'k',clip_on=False)
        ax2.plot([positions[j]],[-0.57],'o',color=colors_r_i_t_m[j],markersize=7,clip_on=False)
    ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=7,clip_on=False)
    ax2.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    plt.gcf().subplots_adjust(bottom=0.3,left=0.12,right=0.99)
    plt.savefig('ffl.pdf',bbox_inches = "tight")
    
    
    
    
def plot_proportion_of_specific_type_ffls_clusters(nr_total_ffl_ffl,nr_specific_ffl_ffl,N):
    import matplotlib
    def arrow_new(self, x, y, dx, dy, **kwargs):
        kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
        kwargs.setdefault('fc', 'black')
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)
        posA = x, y
        posB = x+dx, y+dy
        a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
        self.add_artist(a)
        return a
    
    types_graphic = [list([[0,1],[0,2],[1,2],[3,2],[3,4],[4,2]]),
                     list([[0,1],[0,2],[2,1],[3,2],[3,4],[4,2]]),
                     list([[0,1],[0,2],[1,2],[2,3],[3,4],[2,4]]),
                     list([[0,1],[0,2],[2,1],[3,2],[3,4],[2,4]]),
                     list([[0,1],[2,0],[2,1],[3,2],[3,4],[2,4]]),
                     list([[0,1],[2,0],[2,1],[2,3],[3,4],[2,4]]),
                     list([[0,1],[0,2],[1,2],[3,1],[3,2]]),
                     list([[0,1],[0,2],[1,2],[1,3],[3,2]]),
                     list([[0,1],[0,2],[1,2],[1,3],[2,3]]),
                     list([[1,0],[0,2],[1,2],[1,3],[3,2]]),
                     list([[1,0],[2,0],[1,2],[1,3],[3,2]]),
                     list([[1,0],[2,0],[1,2],[1,3],[2,3]]),
                     list([[0,1],[0,2],[1,2],[2,1]]),
                     list([[0,1],[0,2],[1,2],[1,0]]),
                     list([[0,1],[0,2],[1,2],[1,0],[2,0]])]
    type_graphic = types_graphic[0]
    n_types=15
    order=list(range(n_types))
    all_types = order
    
    cmap = matplotlib.cm.tab20
    
    SHOW_LEGEND = False
    DONT_SHOW_ZERO_FFLS_NETWORKS = True
    LOG=False
    cmap = matplotlib.cm.tab20
    color_ax2 = 'k'
    sorted_sizes = np.array(sorted(zip(nr_total_ffl_ffl,*nr_specific_ffl_ffl/nr_total_ffl_ffl),key=lambda x: (x[0],x[1]))).T
    if DONT_SHOW_ZERO_FFLS_NETWORKS:
        index = list(sorted_sizes[0]>0).index(True)
        sorted_sizes = sorted_sizes[:,index:]
    f,ax=plt.subplots(figsize=(7.5,5) if SHOW_LEGEND else (7.5,3))
    ind = np.arange(sorted_sizes.shape[1])
    bottom = np.zeros(sorted_sizes[order[0]].shape)
    for i in range(n_types):
        ax.bar(ind,sorted_sizes[1+order[i]],bottom=bottom,color=cmap(i),log=LOG,width=0.7)
        bottom +=sorted_sizes[1+order[i]]        
    ax.set_xlim([-1,sorted_sizes.shape[1]])
    ax.set_ylabel('Proportion of specific FFL cluster')
    ax.set_xticks([])
    ax.set_xlabel('Networks sorted by the number of FFL clusters')
    ax2 = ax.twinx()
    ax2.semilogy(ind,sorted_sizes[0],'-',lw=3,color=color_ax2)
    ax2.set_ylabel('Number of FFL clusters',color=color_ax2)
    ax2.tick_params(axis='y', labelcolor=color_ax2)
    delta_y = 0.13
    delta_x = 0.4*((sorted_sizes.shape[1]+1)/7)
    #ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    #ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    #ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
    #ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
    ax.set_ylim([0,1])
    
    if SHOW_LEGEND:
        
        ycenter = 1 + 0.05 + delta_y
        points = [(-1,1),(1,1),(0,0),(-1,-1),(1,-1)] #from top to bottom, left to right
        for j,xcenter in enumerate(np.linspace(-1,sorted_sizes.shape[1],8)[1:-1]):
            head_width = 7
            head_length = 7
            color=cmap(j)
            type_graphic = types_graphic[j]
            for i in range(len(type_graphic)):
                ax.text(xcenter-0.8*delta_x,ycenter,str(all_types[j]),color=color,va='center',ha='center')
                arrow_new(ax,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        
        ycenter = 1 + 0.05 + 3.4*delta_y
        points = [(0,1),(-1,0),(1,0),(0,-1)] #from top to bottom, left to right
        for j,xcenter in enumerate(np.linspace(-1,sorted_sizes.shape[1],8)[1:-1]):
            j=j+6
            head_width = 7
            head_length = 7
            color=cmap(j)
            print(j,color)
            type_graphic = types_graphic[j]
            for i in range(len(type_graphic)):
                ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
                arrow_new(ax,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        
        ycenter = 1 + 0.05 + 4.8*delta_y
        points = [(0,1),(-1,0),(1,0)] #from top to bottom, left to right
        xcenters = np.linspace(-1,sorted_sizes.shape[1],8)[1:-1]
        xcenters = xcenters[0::2]+(xcenters[1]-xcenters[0])/2
        for j,xcenter in enumerate(xcenters):
            j=j+12
            head_width = 7
            head_length = 7
            color=cmap(j)
            print(j,color)
            type_graphic = types_graphic[j]
            for i in range(len(type_graphic)):
                ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
                arrow_new(ax,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
                        
        plt.gcf().subplots_adjust(top=0.57,bottom=0.05)
    plt.savefig('proportion_of_specific_type_ffls_clusters_nice_N%i.pdf' % (N),bbox_inches = "tight")

















if __name__ == '__main__':
    ## Load database
    max_degree = 20
    Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_excluded,models_not_loaded,similar_sets_of_models,n_variables,n_constants,max_degree = load_models_included_in_meta_analysis(max_degree=max_degree,models_to_keep=models_to_keep,models_to_exclude_manually_because_similar_from_same_PID=models_to_exclude_manually_because_similar_from_same_PID)
    N = len(models_loaded)
    
    pmids_models_loaded = [int(''.join([entry for entry in el.split('.')[0].split('_') if (entry[0] in '123456789' and len(entry)>3)])) for el in models_loaded]
    pmids_to_print = list(map(str,pmids_models_loaded))
    dummy = pd.value_counts(pmids_models_loaded)
    PMIDS_with_multiple_included_models  = list(dummy.index)[:sum(dummy>1)]
    for i,PMID in enumerate(pmids_models_loaded):
        if PMID in PMIDS_with_multiple_included_models:
            pmids_to_print[i] += ' (%s)' %  models_loaded[i].split('_')[0].replace('HCC1954 Breast Cell Line ','').replace(' of Drosophila Signaling Pathway','').replace(' of Drosophila Signaling Pathways','').replace(' from the Drosophila Signaling Pathway','').replace(' of Drosophila Signalling Pathways','')
            
    Fs_essential = []
    Is_essential = []
    for i in range(N):
        dummy = can.get_essential_network(Fs[i],Is[i])
        Fs_essential.append(dummy[0])
        Is_essential.append(dummy[1])
    
    ## General analysis
    avg_degrees,avg_essential_degrees = print_basic_summary_and_get_average_degree_per_model(Fs,Is,degrees,degrees_essential,variabless,constantss,N,n_variables,n_constants,max_degree)
    plot_degree_vs_actual_degree(degrees,degrees_essential,N,n_variables,max_degree)
    plot_size_vs_avg_connectivity(n_variables,avg_essential_degrees,N)
    plot_number_of_genes_vs_constants(n_variables,n_constants,N)
    plot_summary_all_networks(n_variables,n_constants,N)
    print_type_of_function(Fs)
    print_common_variables_across_networks(variabless)
    plot_in_and_out_degree(Is,n_variables,degrees,N)
    type_of_each_regulation = compute_type_of_each_regulation_excluding_constants(Fs,n_variables)
    print(pd.value_counts(type_of_each_regulation.type))

    print(pd.value_counts(type_of_each_regulation.type)/type_of_each_regulation.shape[0]*100)
    type_of_each_regulation_by_model = compute_type_of_each_regulation_by_model(type_of_each_regulation,N)
    plt.rcParams.update({'font.size': 9})
    plot_activators_vs_inhibitors_vs_degree(type_of_each_regulation,N,figsize=(2.7,2.4))

    ## Canalization
    depths,layer_structures,all_ncfs,layer_structures_ncf,ks_per_n = get_canalizing_depths(Fs_essential,degrees_essential,n_variables,N)
    write_excel_files_canalizing_depth(ks_per_n,N)

    #new supplementary figure: Figure S3
    import math
    def nchoosek(population, sample):
        "Returns `population` choose `sample`."
        s = max(sample, population - sample)
        assert s <= population
        assert population > -1
        if s == population:
            return 1
        numerator = 1
        denominator = 1
        for i in range(s+1, population + 1):
            numerator *= i
            denominator *= (i - s)
        return numerator/denominator
    def number_canalizing(n):
        return 2*((-1)**n -n - 1) + sum([(-1)**(k+1) * nchoosek(n,k) * 2**(k+1) * 2**(2**(n-k)) for k in range(1,n+1)])    
    upto=10
    log10_proportion_canalizing_expected = []
    for n in range(2,upto+1):
        log10_proportion_canalizing_expected.append(math.log10(number_canalizing(n)) - 2**n * math.log10(2) )
    proportion_canalizing_observed = (np.sum(ks_per_n[:,1:],1)/np.sum(ks_per_n,1))[2:upto+1]
    
    f,ax = plt.subplots(figsize=(3.3,4))
    ax.plot(range(2,upto+1),(proportion_canalizing_observed),'x--',label='observed')
    ax.semilogy(range(2,upto+1),10**np.array(log10_proportion_canalizing_expected),'o:',label='expected')
    ax.legend(loc='best',frameon=False)
    ax.set_xlabel(r'number of inputs, $n$')
    ax.set_ylabel('proportion of canalizing functions\n'+r'(i.e., functions with canalizing depth $k>0$)')
    ax.set_xticks(range(2,upto+1))
    ax.set_yticks([1e-150,1e-120,1e-90,1e-60,1e-30,1])
    ax.set_yticklabels(ax.get_yticklabels()[:-1] + [r'$1$'])
    plt.savefig('prop_canalizing_exp_vs_obs_N%i.pdf' % N,bbox_inches = "tight")    

    def sterling_times_fak_r(n,r):
        return sum([((-1)**i)*nchoosek(r,i)*(r-i)**n for i in range(r+1)])
        
    def sterling_difference(n,r):
        return r**n+sum([(-1)**i*nchoosek(r-1,i-1)*(r-i)**(n-1)*(r**2*1./i-r+n) for i in range(1,r+1)])
    
    def number_ncfs(n,p):
        if n==1 and p==2:
            return (2,0,2)
        if p==2:
            N1=0
        else:
            N1=2**(n-1)*p*(p-2)*(p-1)**n*sum([(p-1)**(r-1)*n*sterling_times_fak_r(n-1,r-1) for r in range(2,n+1)])
        #N2=2**n*p*(p-1)**n*sum([(p-1)**r*(sterling_times_fak_r(n,r)-n*sterling_times_fak_r(n-1,r-1)) for r in range(1,n)])
        N2=2**n*p*(p-1)**n*sum([(p-1)**r*sterling_difference(n,r) for r in range(1,n)])
        return N1+N2,N1,N2

    upto=10
    log10_proportion_nested_canalizing_expected = []
    for n in range(2,upto+1):
        log10_proportion_nested_canalizing_expected.append(math.log10(number_ncfs(n,2)[0]) - 2**n * math.log10(2) )
    proportion_nested_canalizing_observed = (np.diag(ks_per_n)/np.sum(ks_per_n,1))[2:upto+1]
    
    f,ax = plt.subplots(figsize=(3.3,4))
    ax.plot(range(2,upto+1),(proportion_nested_canalizing_observed),'x--',label='observed')
    ax.semilogy(range(2,upto+1),10**np.array(log10_proportion_nested_canalizing_expected),'o:',label='expected')
    ax.legend(loc='best',frameon=False)
    ax.set_xlabel(r'number of inputs, $n$')
    ax.set_ylabel('proportion of nested canalizing functions\n'+r'(i.e., functions with canalizing depth $k=n$)')
    ax.set_xticks(range(2,upto+1))
    ax.set_yticks([1e-300,1e-240,1e-180,1e-120,1e-60,1])
    ax.set_yticklabels(ax.get_yticklabels()[:-1] + [r'$1$'])
    plt.savefig('prop_nested_canalizing_exp_vs_obs_N%i.pdf' % N,bbox_inches = "tight")    

    res_from_get_number_of_layers = get_number_of_layers(Fs_essential)
    for degree in [3,4,5]:
        get_canalizing_input_and_output_distribution_of_NCFs(res_from_get_number_of_layers,N,degree=degree)
    get_count_of_specific_NCF_layer_structures(layer_structures_ncf,N)
    ss_per_n = get_symmetry_groups_of_networks(Fs,degrees,degrees_essential,N,upto=min(10,max_degree))
    ss_per_n_random,ss_per_n_random_canalizing,ss_per_n_random_NCF,ss_per_n_imputed = get_random_expectations_symmetry_groups(ks_per_n)
    plot_symmetry_distributions(ss_per_n,ss_per_n_imputed,ss_per_n_random,N)

    plot_canalizing_strength_of_noncanalizing_functions(N,max_degree=6,nsim=1000)
    plot_effectiveness_of_noncanalizing_functions(N,max_degree=6,nsim=1000)

    ## FBLs
    max_loop=6
    nr_loops,nr_pos_loops,nr_neg_loops,nr_unknown_loops,nr_notreal_loops,nr_specific_k_loops,nr_real_loops,all_types,all_loops = compute_all_FBLs(Fs,Is,degrees,constantss,max_degree,max_loop=max_loop)
    plot_number_of_pos_neg_and_unknown_FBLs(nr_real_loops,nr_pos_loops,nr_neg_loops,nr_unknown_loops,N,max_loop=max_loop)
    plot_total_count_of_FBLs(nr_pos_loops,nr_neg_loops,nr_unknown_loops,N,max_loop=max_loop)
    plot_total_count_of_specific_FBLs(nr_specific_k_loops,nr_unknown_loops,N,max_loop=max_loop,LEGEND=False)
    
    prop_pos_separate,prop_pos_separate_by_scc,prop_pos_specific_scc_of_loop,prop_pos_specific_scc_of_loop2 = get_null_expectations_fbl(Fs,Is,variabless,all_types,all_loops)
    plot_total_count_of_pos_vs_neg_FBLs(nr_specific_k_loops,nr_unknown_loops,prop_pos_separate,prop_pos_specific_scc_of_loop2,N,max_loop=max_loop)
    plot_pos_vs_neg_regulations_in_FBLs_best(Fs,nr_specific_k_loops,all_loops,prop_pos_separate,prop_pos_specific_scc_of_loop,N)
    
    
    ## FFLs
    nr_ffls,nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls,nr_specific_ffls,nr_real_ffls,nr_notreal_ffls = compute_all_FFLs(Fs,Is,degrees,constantss,max_degree,N)
    plot_number_of_coherent_vs_incoherent_ffls_by_model(nr_real_ffls,nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls,N)
    plot_proportion_of_coherent_vs_incoherent_ffls_by_model(nr_ffls,nr_real_ffls,nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls,N)
    expected_number_per_type = get_null_expectations_ffl(nr_specific_ffls,prop_pos_separate,N)
    order_ffls = plot_total_count_of_specific_ffls(nr_specific_ffls,nr_unknown_ffls,expected_number_per_type,N)
    plot_proportion_of_specific_ffls_per_model(nr_real_ffls,nr_specific_ffls,nr_unknown_ffls,N)
    
    #compute_layer_regulator_vs_layer_intermediate_in_ffls(Fs,Is,constantss)
    res,_ = compute_effectiveness_regulator_vs_effectiveness_intermediate_in_ffls(Fs,Is,constantss,upto=min(7,max_degree))
    plot_effectiveness_regulator_vs_effectiveness_intermediate_in_ffls(res,N,upto=min(7,max_degree),figsize=(4.2,3))
    

    ## FFL clusters
    nr_total_ffl_ffl,nr_specific_ffl_ffl = compute_all_FFL_clusters(Is,constantss,n_variables,degrees,max_degree,N)
    print(np.sum(nr_specific_ffl_ffl),np.sum(nr_specific_ffl_ffl[:6]),np.sum(nr_specific_ffl_ffl[6:12]),np.sum(nr_specific_ffl_ffl[12:]))
    print(np.sum(nr_specific_ffl_ffl[:6])/np.sum(nr_specific_ffl_ffl)*100,np.sum(nr_specific_ffl_ffl[6:12])/np.sum(nr_specific_ffl_ffl)*100,np.sum(nr_specific_ffl_ffl[12:])/np.sum(nr_specific_ffl_ffl)*100)
    plt.rcParams.update({'font.size': 12})
    plot_total_count_of_specific_type_ffls_clusters(nr_total_ffl_ffl,nr_specific_ffl_ffl,N)
    plot_proportion_of_specific_type_ffls_clusters(nr_total_ffl_ffl,nr_specific_ffl_ffl,N)
    plot_total_count_of_specific_type_ffls_clusters_firsttwelve(nr_total_ffl_ffl,nr_specific_ffl_ffl,N)
    plot_total_count_of_specific_type_ffls_clusters_lastthree(nr_total_ffl_ffl,nr_specific_ffl_ffl,N)
    plt.rcParams.update({'font.size': 9})
    
    


    
    ## Dynamics - derrida values
    derrida_values = get_derrida_values(Fs,Is,degrees,max_degree,nsim=10000)
    plot_derrida_values_vs_network_size(derrida_values,n_variables,max_degree,N)
    plot_derrida_values_vs_network_size(derrida_values,avg_essential_degrees,max_degree,N)
    
    
    ## write everything in one output file
    general = np.c_[pmids_to_print,models_loaded,n_variables,n_constants,avg_degrees,avg_essential_degrees,type_of_each_regulation_by_model]
    general_labels = (['filename','models_loaded','number of variables','number of constants','average in-degree','average essential in-degree'] + ['number of regulations that are']*4)
    general_labels2 = ([''] * 6 + ['positive','negative','conditional','not essential'])
    A = pd.DataFrame(np.r_[[general_labels2],general],columns = general_labels)
    A.to_csv('general_info_N%i.csv' % N)
    
    
    general = np.c_[pmids_to_print,
                    nr_ffls,nr_coherent_ffls,nr_incoherent_ffls,nr_unknown_ffls,nr_real_ffls,nr_notreal_ffls,nr_specific_ffls[order_ffls,:].T]
    general_labels = (['filename'] +
                      ['number of feed-forward loops'] * 14)
    general_labels2 = ([''] + 
                       ['total','coherent','incoherent','unknown type (i.e., including a conditional edge)','real (i.e., all three edges essential)','not real (including a non-essential edge)'] + 
                       ['type '+str(i+1) for i in range(8)])
    A = pd.DataFrame(np.r_[[general_labels2],general],columns = general_labels)
    A.to_csv('ffl_info_N%i.csv' % N)
    
    
    general = np.c_[pmids_to_print]
    general_labels = ['filename']
    general_labels2 = ['']
    for i in range(max_loop):
        general = np.c_[general,nr_loops[i],nr_pos_loops[i],nr_neg_loops[i],nr_unknown_loops[i],nr_real_loops[i],nr_notreal_loops[i]]
        general_labels.extend(['number of '+str(i+1)+'-loops'] * (6+(i+2)))
        general_labels2.extend(['total','positive','negative','unknown type (i.e., including a conditional edge)','real (i.e., all edges essential)','not real (including a non-essential edge)'])
        for j in range(i+2):
            general = np.c_[general,nr_specific_k_loops[i,j]]
            general_labels2.append('%i+ %i-' % (i+1-j,j))
    A = pd.DataFrame(np.r_[[general_labels2],general],columns = general_labels)
    A.to_csv('fbl_info_N%i.csv' % N)
    
    
    
    general = np.c_[pmids_to_print,
                    nr_total_ffl_ffl,nr_specific_ffl_ffl.T]
    general_labels = (['filename'] +
                      ['number of feed-forward loop clusters'] * 16)
    general_labels2 = (['','total'] + 
                       ['type '+str(i+1) for i in range(15)])
    A = pd.DataFrame(np.r_[[general_labels2],general],columns = general_labels)
    A.to_csv('ffl_cluster_info_N%i.csv' % N)
    
    
    def flatten(l):
        return [item for sublist in l for item in sublist]
    
    
    #load kingdom information to stratify by kingdom
    general_info_csv = pd.read_excel('general_info_N122_models.xlsx')   
    pmids_in_excel = np.array(general_info_csv['PMID (model)'][1:])
    kingdoms_in_excel = np.array(general_info_csv['Organism'][1:])
    for i in range(N):
        if ':' in kingdoms_in_excel[i]:
            kingdoms_in_excel[i] = kingdoms_in_excel[i].split(':')[0]
    pmids_to_print_dict = dict(zip(pmids_to_print,list(range(N))))
    kingdoms = ['']*N
    for i in range(N):
        kingdoms[pmids_to_print_dict[str(pmids_in_excel[i])]] = kingdoms_in_excel[i]
    kingdoms = np.array(kingdoms)
    print(pd.value_counts(kingdoms))
    kingdoms_of_interest = ['Animal','Bacteria','Fungi','Plant']
    
    for kingdom_of_interest in kingdoms_of_interest:
        print(kingdom_of_interest)
        print()
        which = kingdoms == kingdom_of_interest
        avg_degrees,avg_essential_degrees = print_basic_summary_and_get_average_degree_per_model(list(np.array(Fs)[which]),list(np.array(Is)[which]),list(np.array(degrees)[which]),list(np.array(degrees_essential)[which]),list(np.array(variabless)[which]),list(np.array(constantss)[which]),sum(which),n_variables[which],n_constants[which],max_degree)
        print()
        print()

    sum_types_of_regulation_by_kingdom = []
    for kingdom_of_interest in kingdoms_of_interest:

        print(kingdom_of_interest)
        print()
        which = kingdoms == kingdom_of_interest
        type_of_each_regulation = compute_type_of_each_regulation_excluding_constants(list(np.array(Fs)[which]),n_variables[which])
        plot_activators_vs_inhibitors_vs_degree(type_of_each_regulation,kingdom_of_interest,figsize=(4,2.4))
        print()
        print()

        sum_types_of_regulation_by_kingdom.append([])
        #for label in ['increasing','decreasing','not monotonic','not essential']:
        for label in ['increasing','decreasing','not monotonic']:
            sum_types_of_regulation_by_kingdom[-1].append(sum(type_of_each_regulation.type==label))
    
    sum_types_of_regulation_by_kingdom = np.array(sum_types_of_regulation_by_kingdom)
    
    f,ax = plt.subplots()
    for i in range(3):
        ax.bar(range(4),sum_types_of_regulation_by_kingdom[:,i]/np.sum(sum_types_of_regulation_by_kingdom,1),bottom = np.sum(sum_types_of_regulation_by_kingdom[:,:i],1)/np.sum(sum_types_of_regulation_by_kingdom,1))
    
        

    
    max_loop=6
    nr_loops,nr_pos_loops,nr_neg_loops,nr_unknown_loops,nr_notreal_loops,nr_specific_k_loops,nr_real_loops,all_types,all_loops = compute_all_FBLs(Fs,Is,degrees,constantss,max_degree,max_loop=max_loop)
    prop_pos_separate,prop_pos_separate_by_scc,prop_pos_specific_scc_of_loop,prop_pos_specific_scc_of_loop2 = get_null_expectations_fbl(Fs,Is,variabless,all_types,all_loops)

    for kingdom_of_interest in kingdoms_of_interest:
        which = kingdoms == kingdom_of_interest
        expected_number_per_type = get_null_expectations_ffl(nr_specific_ffls[:,which],prop_pos_separate[which],sum(which))
        order_ffls = plot_total_count_of_specific_ffls(nr_specific_ffls[:,which],nr_unknown_ffls[which],expected_number_per_type,kingdom_of_interest)

    for kingdom_of_interest in kingdoms_of_interest:
        which = kingdoms == kingdom_of_interest    
        prop_pos_separate,prop_pos_separate_by_scc,prop_pos_specific_scc_of_loop,prop_pos_specific_scc_of_loop2 = get_null_expectations_fbl(list(np.array(Fs)[which]),list(np.array(Is)[which]),list(np.array(variabless)[which]),list(np.array(all_types)[which]),list(np.array(all_loops)[which]))
        plot_total_count_of_pos_vs_neg_FBLs(nr_specific_k_loops[:,:,which],nr_unknown_loops[:,which],prop_pos_separate,prop_pos_specific_scc_of_loop2,kingdom_of_interest,max_loop=max_loop)
        plot_pos_vs_neg_regulations_in_FBLs_best(list(np.array(Fs)[which]),nr_specific_k_loops[:,:,which],all_loops,prop_pos_separate,prop_pos_specific_scc_of_loop,kingdom_of_interest)
        
    
    
    variabless_simple = [[el.lower().replace('_','').replace('.','').replace('kappa','k') for el in variables] for variables in variabless]
    
    
    all_variables = []
    for i in range(len(variabless)):
        #all_variables.extend(list(map(str.lower,variabless[i])))
        all_variables.extend(variabless_simple[i])
    all_variables = np.array(all_variables)  
    
    
    kingdoms_of_interest = ['Animal','Bacteria','Fungi','Plant','multiple kingdoms']

    all_unique_variables = list(set(all_variables))
    all_unique_variables_dict = dict(zip(all_unique_variables,range(len(all_unique_variables))))
    count = np.zeros((len(all_unique_variables_dict),len(kingdoms_of_interest)),dtype=int)
    for j,kingdom_of_interest in enumerate(kingdoms_of_interest):
        which = kingdoms == kingdom_of_interest 
        for i in range(len(variabless_simple)):
            if which[i]:
                for v in variabless_simple[i]:
                    count[all_unique_variables_dict[v],j] += 1
    A = pd.DataFrame(np.c_[all_unique_variables,np.sum(count,1),count])
    A.to_excel('all_variables_with_count.xlsx')
    
    
        
        
            
        
    
    
    
    
    
    
    
    
    #additional figures on criticality and NCFs
    IR = [] # IR = input redundancy (Gates et al., PNAS, 2021)
    degs = []
    for i,F in enumerate(Fs):
        IR.append([])
        degs.append([])
        for f in F[:n_variables[i]]:
            degs[-1].append(0 if len(f)==0 else int(np.log2(len(f))))
            if len(f)>0 and len(f)<=2**8:
                IR[-1].append(can.get_input_redundancy(f))
            else:
                IR[-1].append(np.nan)
                
    EC = [] # EC = effective connectivity (Gates et al., PNAS, 2021)
    for F in Fs:
        EC.append([])
        for f in F:
            if len(f) == 0:
                n=0
            else:
                n = int(np.log2(len(f)))
            if n>0 and n<=12:
                EC[-1].append( sum(can.get_edge_effectiveness(f,)) )
            else:
                EC[-1].append( np.nan )

    hamming_weight = [[sum(f) for f in F] for F in Fs_essential]
    deg_essential = [[0 if len(f)==0 else int(np.log2(len(f))) for f in F] for F in Fs_essential]
    bias = [[hw/2**d*(1-hw/2**d) for hw,d in zip(hamming_weight[i],degrees[i])] for i in range(N)]
    
    #IR_flat = np.array(flatten(IR))
    hamming_weight_flat = np.array(flatten(hamming_weight))
    degrees_essential_flat = np.array(flatten(degrees_essential))
    degs_flat = np.array(flatten(degs))
    degrees_flat = np.array(flatten(degrees))
    bias_flat = np.array(flatten(bias))
    
    
    
    plt.plot(flatten(degrees_essential),flatten(bias),'o')
    
    prop_ncf = np.array([np.mean([depths[i][j]==degrees_essential[i][j] for j in range(n_variables[i])]) for i in range(N)])
    mean_bias = np.array([np.nanmean(el) for el in bias])
    mean_EC = np.array([np.nanmean(el) for el in EC])
    mean_IR = np.array([np.nanmean(el) for el in IR])
    mean_deg_essential = [np.mean(degrees_essential[i][:n_variables[i]]) for i in range(N)]
    derrida = np.array(derrida_values)
    
    f,ax = plt.subplots()
    for VAL in [0,1]:
        ax.plot([np.mean(degrees_essential[i][:n_variables[i]]) for i in range(N) if all_ncfs[i]==VAL],[derrida_values[i] for i in range(N) if all_ncfs[i]==VAL],'o')

    f,ax = plt.subplots()
    ax.scatter(mean_deg_essential,derrida_values,c=prop_ncf)

    f,ax = plt.subplots()
    ax.scatter(mean_deg_essential,derrida_values,c=mean_bias)

    f,ax = plt.subplots()
    ax.scatter(mean_deg_essential,derrida_values,c=mean_IR)

    f,ax = plt.subplots()
    ax.scatter(mean_deg_essential,mean_bias,c=derrida_values)

    f,ax = plt.subplots()
    ax.scatter(mean_deg_essential,prop_ncf,c=derrida_values)


    degs_flat
    depths_flat = np.array(flatten(depths))
    layer_structures_flat = flatten(layer_structures)

    n=3
    indices = np.bitwise_and(degs_flat==n,depths_flat==n)
    pd.value_counts(np.array(layer_structures_flat,dtype='object')[indices])
    
    nsim=1000
    res = []
    for _ in range(nsim):
        f = can.random_k_canalizing(n, n)
        res.append(can.get_layerstructure_given_canalizing_outputs_and_corefunction(can.find_layers(f)[3],[1]))
    res = np.array(res)
    pd.value_counts(res)
    

    def activities_NCF(layerstructure_ncf):
        '''Theorem 3.7'''
        r = len(layerstructure_ncf)
        n = sum(layerstructure_ncf)
        phi = np.zeros(r+2)
        psi = np.zeros(r+1)
        psi[r::-2] = 1
        for i in range(r-1,0,-1):
            phi[i] = phi[i+2] + sum([0.5**(sum(layerstructure_ncf[:i])+s) for s in range(0,layerstructure_ncf[i])])
        return phi[1:-1] + 1/2**(n-1) * psi[1:]
    
    def average_sensitivity_NCF(layerstructure_ncf):
        act = activities_NCF(layerstructure_ncf)
        return np.dot(act,layerstructure_ncf)
    
    import matplotlib
    cmap = matplotlib.cm.Set1
    colors = [cmap(i) for i in range(6)]
    markers = ['P','X','o','v','s']
    f,ax = plt.subplots(figsize=(3.7,2.5))
    means=[]
    for i,n in enumerate([3,4,5,6]):
        indices = np.bitwise_and(degs_flat==n,depths_flat==n)
        n_NCFs = sum(indices)
        dummy = pd.value_counts(np.array(layer_structures_flat,dtype='object')[indices])
        props = []
        derridas = []
        for el,val in zip(dummy.index,dummy.values):
            derridas.append( average_sensitivity_NCF(el) )
            props.append( val / n_NCFs)
        ax.plot(derridas,props,marker=markers[i],color=colors[i],ls='none')
        mean = np.dot(props,derridas)
        means.append(mean)
    [y1,y2] = ax.get_ylim()
    for i,n in enumerate([3,4,5,6]):    
        ax.plot([means[i],means[i]],[-10,10],marker=markers[i],color=colors[i],label=str(n))
    ax.set_ylim(y1,y2)
    ax.legend(loc='best',frameon=False,ncol=2,title='k')
    ax.set_xlabel('Average sensitivity of NCF')
    ax.set_ylabel('Relative proportion (per k)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig('prop_different_NCFs_N%i.pdf' % N,bbox_inches = "tight")    
            
    #random exp.
    n=5
    avgs = []
    derridas = []
    for w in range(1,2**(n-1),2):
        ls = can.get_layerstructure_of_an_NCF_given_its_Hamming_weight(n,w)[1]
        avgs.append(average_sensitivity_NCF(ls))
        derridas.append(np.dot(ls,avgs[-1]))
        

    f,ax = plt.subplots(figsize=(2.8,2.5))
    im = ax.scatter(mean_deg_essential,derrida_values,c=prop_ncf,s=15)
    cbar = plt.colorbar(im)
    cbar.set_label('Proportion of NCFs')    
    ax.set_xlabel('Average essential degree')
    ax.set_ylabel('Mean average sensitivity')
    plt.savefig('derrida_vs_ess_degree_color_by_prop_ncf_N%i.pdf' % N,bbox_inches = "tight")    
    
    A = pd.DataFrame(np.c_[mean_deg_essential,derrida_values,prop_ncf],columns = ['mean deg essential','derrida','prop NCF'])
    A.to_excel('out2.xlsx')
    
    import scipy.stats as stats
    stats.pearsonr(np.array(mean_deg_essential)[np.isnan(derrida_values)==0],np.array(derrida_values)[np.isnan(derrida_values)==0])
    stats.pearsonr(np.array(n_variables)[np.isnan(derrida_values)==0],np.array(derrida_values)[np.isnan(derrida_values)==0])
    stats.pearsonr(np.array(prop_ncf)[np.isnan(derrida_values)==0],np.array(derrida_values)[np.isnan(derrida_values)==0])
    
        
    import matplotlib
    
    cov = [np.cov(k,p)[0,1] for k,p in zip(deg_essential,bias)]
    correl = [stats.pearsonr(k,p)[0] for k,p in zip(deg_essential,bias)]
    PEARSON = True
    func = stats.pearsonr if PEARSON else stats.spearmanr
    notnan = np.isnan(derrida_values)==False
    ys =[n_variables,prop_ncf,mean_bias,mean_deg_essential,mean_EC,cov,
         np.array(mean_bias)*np.array(mean_deg_essential),
         np.array(mean_bias)*np.array(mean_deg_essential)+cov,
         np.array(mean_bias)*np.array(mean_EC),
         derrida_values]
    labels = 'network size,proportion NCF,<p(1-p)>,<$K$>,<$K_e$>,Cov(p(1-p);K),<K><p(1-p)>,<K><p(1-p)>+Cov,<Ke><p(1-p)>,mean average sensitivity'.replace('_','\n').split(',')
    labels[2] = r'<$p(1-p)$>'
    labels[3] = r'<$K$>'
    labels[4] = r'<$K_e$>'
    labels[5] = r'Cov($p(1-p),K$)'
    labels[-3] = r'<$K$><$p(1-p)$> + Cov'
    labels[-2] = r'<$K_e$><$p(1-p)$>'
    n_ys = len(ys)
    correl_mat = np.ones((n_ys,n_ys))
    for i in range(n_ys):
        for j in range(i+1,n_ys):
            correl_mat[i,j] = func(np.array(ys[i])[notnan],np.array(ys[j])[notnan])[0]
            correl_mat[j,i] = correl_mat[i,j] 
    f,ax = plt.subplots(figsize=(3.4,2.5))   
    im = ax.imshow(correl_mat,cmap=matplotlib.cm.RdBu,vmin=-1,vmax=1)
    cbar = plt.colorbar(im)
    cbar.set_label(('Pearson' if PEARSON else 'Spearman') + ' correlation')
    ax.set_xticks(range(n_ys))
    ax.set_xticklabels(labels,rotation=90)    
    ax.set_yticks(range(n_ys))
    ax.set_yticklabels(labels)     
    plt.savefig(('Pearson' if PEARSON else 'Spearman')+'_derrida_N%i.pdf' % N,bbox_inches = "tight")    
    
    
    n=3
    indices = np.bitwise_and(degs_flat==n,depths_flat==n)
    pd.value_counts(np.array(layer_structures_flat,dtype='object')[indices])
    
    res = []
    n_obs_per_n = {}
    for n in range(3,8):

        indices = np.bitwise_and(degs_flat==n,depths_flat==n)
        dummy = pd.value_counts(np.array(layer_structures_flat,dtype='object')[indices])
        dummy_index = list(dummy.index)
        dummy_values = list(dummy.values)
        sum_dummy_values = sum(dummy_values)
        n_obs_per_n.update({n:sum_dummy_values})
        
        for w in range(1,2**(n-1),2):
            r,ls = can.get_layerstructure_of_an_NCF_given_its_Hamming_weight(n,w)
            f = can.random_k_canalizing_with_specific_layerstructure(n,ls)
            effective_connectivity = sum(can.get_edge_effectiveness(f,n))
            average_sensitivity = average_sensitivity_NCF(ls)
            bias_1_m_bias = w/2**n * (1 - w/2**n)
            
            try:
                n_obs = dummy_values[dummy_index.index(ls)]
            except ValueError:
                n_obs = 0
            prop_obs = n_obs/sum_dummy_values
            
            res.append([n,n_obs,100*prop_obs,w,effective_connectivity,bias_1_m_bias,average_sensitivity,r,', '.join(list(map(str,ls)))])
    
    A = pd.DataFrame(res,columns = 'n,n_obs,n_obs (in %),w,effective_connectivity,bias_1_m_bias,average_sensitivity,number of layers,layer structure'.split(','))
    A.to_excel('out_ncf.xlsx')
    
    import scipy
    width=0.7
    f,ax = plt.subplots()
    colors = [matplotlib.cm.Set1(i) for i in range(8)]
    for n in range(3,8):
        props = []
        exp_props = np.array([scipy.special.binom(n-2,k) for k in range(n-1)])/2**(n-2)
        for r in range(1,n):
            which = np.bitwise_and(A.n==n,A['number of layers']==r)
            props.append(sum(A.n_obs[which])/n_obs_per_n[n])
            ax.bar([2*n+1],props[-1],bottom = sum(props[:r-1]),color = colors[r-1])
            ax.bar([2*n],exp_props[r-1],bottom = sum(exp_props[:r-1]),color = colors[r-1])
        



    f,ax = plt.subplots(figsize=(4.2,1.8))
    width=0.7
    epsilon = 0.1
    colors = [matplotlib.cm.Set2(i) for i in range(8)]
    max_n = 7
    
    for i,n in enumerate(range(3,max_n+1)):
        props = []
        exp_props = np.array([scipy.special.binom(n-2,k) for k in range(n-1)])/2**(n-2)
        for r in range(1,n):
            which = np.bitwise_and(A.n==n,A['number of layers']==r)
            props.append(sum(A.n_obs[which])/n_obs_per_n[n])
            ax.bar([2*i+width/2+epsilon],props[-1],bottom = sum(props[:r-1]),color = colors[r-1],width=width)
            ax.bar([2*i-width/2-epsilon],exp_props[r-1],bottom = sum(exp_props[:r-1]),color = colors[r-1],width=width)

    lw=1
    for i,n in enumerate(range(3,max_n+1)):
        for y in [0,1]:
            ax.plot([2*i-width-epsilon,2*i-width-epsilon+width],[y,y],'k--',lw=lw)
            ax.plot([2*i+epsilon,2*i+epsilon+width],[y,y],'k-',lw=lw)
        for dx in [0,width]:
            ax.plot([2*i-width-epsilon+dx,2*i-width-epsilon+dx],[0,1],'k--',lw=lw)
            ax.plot([2*i+epsilon+dx,2*i+epsilon+dx],[0,1],'k-',lw=lw)

    x1,x2 = ax.get_xlim()
    
    height_rectangle = 0.07
    width_rectangle = 0.5
    start_x = -0.4
    text_offset = 0.3
    skip_x = 2
    ax.text(start_x+width_rectangle+text_offset+skip_x*(1%3)-0.25,1.32,'number of layers of NCF',ha='center',va='center')
    for j in range(max_n-1):
        ax.add_patch(matplotlib.patches.Rectangle([start_x+skip_x*(j%3),1.18 if j//3==0 else 1.05],width_rectangle,height_rectangle,color=colors[j],clip_on=False))
        ax.text(start_x+width_rectangle+text_offset+skip_x*(j%3),(1.18 if j//3==0 else 1.05)+height_rectangle/2.5,str(j+1),ha='left',va='center',clip_on=False)

    start_x = 5.7
    
    for j,(ls,label) in enumerate(zip(['--','-'],['expected','observed'])):
        ax.text(start_x+width_rectangle+text_offset,1.18 - 0.13*j+height_rectangle/2.5,label,ha='left',va='center',clip_on=False)
        ax.plot([start_x,start_x+width_rectangle],[1.18 - 0.13*j,1.18 - 0.13*j],'k',ls=ls,lw=1,clip_on=False)
        ax.plot([start_x,start_x+width_rectangle],np.array([1.18 - 0.13*j,1.18 - 0.13*j])+height_rectangle,'k',ls=ls,lw=1,clip_on=False)
        ax.plot([start_x,start_x],np.array([1.18 - 0.13*j,1.18 - 0.13*j+height_rectangle]),'k',ls=ls,lw=1,clip_on=False)
        ax.plot([start_x+width_rectangle,start_x+width_rectangle],np.array([1.18 - 0.13*j,1.18 - 0.13*j+height_rectangle]),'k',ls=ls,lw=1,clip_on=False)
    
    ax.set_xlim([x1,x2])    
    ax.set_ylim([0,1])
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('none')
    ax.xaxis.set_ticks(2*np.arange(len(range(3,8))))
    ax.xaxis.set_ticklabels(list(map(str,range(3,8))))
    ax.set_ylabel('Proportion')
    ax.set_xlabel('In-degree')
    plt.savefig('types_of_NCFs_N%i.pdf' % N,bbox_inches = "tight")    




