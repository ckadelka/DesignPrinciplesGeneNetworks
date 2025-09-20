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




if __name__ == '__main__':
    ## Load database
    max_degree = 20 # load all rules with 20 or less regulators
    Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_excluded,models_not_loaded,similar_sets_of_models,n_variables,n_constants,max_degree = load_models_included_in_meta_analysis(max_degree=max_degree,models_to_keep=models_to_keep,models_to_exclude_manually_because_similar_from_same_PID=models_to_exclude_manually_because_similar_from_same_PID)
    N = len(models_loaded)
    
            
    Fs_essential = []
    Is_essential = []
    for i in range(N):
        dummy = can.get_essential_network(Fs[i],Is[i])
        Fs_essential.append(dummy[0])
        Is_essential.append(dummy[1])
    avg_degrees,avg_essential_degrees = print_basic_summary_and_get_average_degree_per_model(Fs,Is,degrees,degrees_essential,variabless,constantss,N,n_variables,n_constants,max_degree)
    pmids_models_loaded = [int(''.join([entry for entry in el.split('.')[0].split('_') if (entry[0] in '123456789' and len(entry)>3)])) for el in models_loaded]
    pmids_to_print = list(map(str,pmids_models_loaded))
    dummy = pd.value_counts(pmids_models_loaded)
    PMIDS_with_multiple_included_models  = list(dummy.index)[:sum(dummy>1)]
    for i,PMID in enumerate(pmids_models_loaded):
        if PMID in PMIDS_with_multiple_included_models:
            pmids_to_print[i] += ' (%s)' %  models_loaded[i].split('_')[0].replace('HCC1954 Breast Cell Line ','').replace(' of Drosophila Signaling Pathway','').replace(' of Drosophila Signaling Pathways','').replace(' from the Drosophila Signaling Pathway','').replace(' of Drosophila Signalling Pathways','')
                
    
    
    ## Compute strongly connected components (SCCs)
    res = []
    sccs = []
    n_sccs = []
    len_sccs = []

    import networkx as nx
    all_dict_sccs = []
    for i,(F,I,var) in enumerate(zip(Fs_essential,Is_essential,variabless)): 
        n_var = len(var)
        
        G = can.generate_networkx_graph_from_edges(I,len(I))
        sccs.append([list(scc) for scc in nx.strongly_connected_components(G)])
        n_sccs.append( len(sccs[-1]) )
        len_sccs.append( np.sort(list(map(len,sccs[-1])))[::-1])
        dict_sccs = {}
        for i in range(n_sccs[-1]):
            for el in sccs[-1][i]:
                dict_sccs.update({el:i})
        all_dict_sccs.append(dict_sccs)
        
        
        
    ## Compute number of non-trivial SCCs (i.e., those containing more than one node)
    counter=0
    counter_no_loop=0
    indices = []
    indices_no_loop = []
    n_nontrivial_sccs = []
    for i in range(N):
        if len(len_sccs[i])>1 and len_sccs[i][1]>1:
            counter+=1
            indices.append(i)
            print(i,len(len_sccs[i][len_sccs[i]>1]),len_sccs[i][len_sccs[i]>1])

        if len_sccs[i][0]==1:
            counter_no_loop+=1
            indices_no_loop.append(i)
        n_nontrivial_sccs.append(len(len_sccs[i][len_sccs[i]>1]))
    n_nontrivial_sccs = np.array(n_nontrivial_sccs)
    
    print(sum(n_nontrivial_sccs==0))
    print(sum(n_nontrivial_sccs>1))
    
    
    ## Spearman correlations
    import scipy.stats as stats
    print(stats.spearmanr(n_variables,n_nontrivial_sccs))
    print(stats.spearmanr(avg_essential_degrees,n_nontrivial_sccs))
    
    ## Logistic regression on the binary variable: model has more than on non-trivial SCC
    import statsmodels.api as sm
    X = np.c_[n_variables,avg_essential_degrees,np.ones(N)] 
    y = n_nontrivial_sccs>1
    logit_model=sm.Logit(y,X)   
    result=logit_model.fit()
    print(result.summary2())
    

    
    
    
    
    ## Plot the directed acyclic graph (DAG)
    n_generationss = []
    n_nontrivial_sccs_indices = n_nontrivial_sccs[indices]
    n_sccs = [len(sccs[i]) for i in indices]
    for i in indices:
        scc_dict = {}
        for j,s in enumerate(sccs[i]):
            for el in s:
                scc_dict.update({el:j})

                
        outdegrees = [0 for j in range(len(Fs[i]))]
        for regulators in Is_essential[i][:n_variables[i]]:
            for regulator in regulators:
                try:
                    outdegrees[regulator] += 1
                except IndexError: #if regulator is a constant
                    pass
        
        len_sccs = np.array(list(map(len,sccs[i])))
        types_sccs = []
        for ii,scc in enumerate(sccs[i]):
            if len_sccs[ii]>1:
                types_sccs.append(2)
            elif scc[0] >= n_variables[i]: #input
                types_sccs.append(-1)
            elif outdegrees[scc[0]] == 0: #output
                types_sccs.append(1)
            else: #transient node
                types_sccs.append(0)
        types_sccs = np.array(types_sccs)
            
                
        
        dag = set()
        for ii,I in enumerate(Is_essential[i][:n_variables[i]]):
            for regulator in I:
                edge = (scc_dict[regulator],scc_dict[ii])
                if len(sccs[i])>40:
                    if edge[0]!=edge[1] and regulator<n_variables[i] and outdegrees[ii]>0:
                        dag.add(edge)
                else:
                    if edge[0]!=edge[1]:
                        dag.add(edge)                    
        nodes_in_dag = list(set.union(*[set(el) for el in dag]))
        
        labels = dict(zip(nodes_in_dag,[str(len_sccs[node]) if len_sccs[node]>1 else '' for node in nodes_in_dag]))
        g = nx.DiGraph(dag)
        top_to_bottom = [sorted(generation) for generation in nx.topological_generations(g)]
        pos = dict()
        n_generations = 0
        for ii,generation in enumerate(nx.topological_generations(g)):
            len_generation = len(generation)
            if len_generation==1:
                pos.update( {generation[0] : (0,-ii)})
            else:
                delta_x = min(len_generation,7)/7
                for jj,node in enumerate(generation):
                    #pos.update( {node : (-0.5+jj/(len_generation-1),-ii)})
                    pos.update( {node : (-delta_x/2+delta_x*jj/(len_generation-1),-ii)})
            n_generations+=1
        n_generationss.append(n_generations)
        f,ax = plt.subplots()
        colors = np.array(['#555555','#ffaaaa','#555555','#ff7777'])[types_sccs[np.array(g.nodes)]+1]
        nx.draw_networkx(g,labels=labels,pos=pos,node_color = colors,#edgecolors='#878787',
                         node_size=200*np.log(1.3*len_sccs[np.array(g.nodes)]))
        #ax.set_title(str(i)+' : '+models_loaded[i])
        ax.set_title(pmids_to_print[i])
        plt.axis('off')
        plt.savefig('sccs_%s.pdf' % (pmids_to_print[i].replace(' ','_').replace('-','_').replace('(','_').replace(')','')),bbox_inches = "tight")    

    sorted_ = np.array(sorted(range(len(indices)), key = lambda x: (n_generationss[x], n_sccs[x])))
    n_generationss = np.array(n_generationss)[sorted_]
    indices_sorted = np.array(indices)[sorted_]
    
    
    ## Generate one big plot of the DAG structure all 30 models with more than one non-trivial component
    f,ax = plt.subplots(nrows=6,ncols=5,figsize=(12,16),gridspec_kw={'height_ratios': n_generationss.reshape((6,5)).max(1)})
    ax = ax.flatten()
    for iii,i in enumerate(indices_sorted):
        scc_dict = {}
        for j,s in enumerate(sccs[i]):
            for el in s:
                scc_dict.update({el:j})

                
        outdegrees = [0 for j in range(len(Fs[i]))]
        for regulators in Is_essential[i][:n_variables[i]]:
            for regulator in regulators:
                try:
                    outdegrees[regulator] += 1
                except IndexError: #if regulator is a constant
                    pass
        
        len_sccs = np.array(list(map(len,sccs[i])))
        types_sccs = []
        for ii,scc in enumerate(sccs[i]):
            if len_sccs[ii]>1:
                types_sccs.append(2)
            elif scc[0] >= n_variables[i]: #input
                types_sccs.append(-1)
            elif outdegrees[scc[0]] == 0: #output
                types_sccs.append(1)
            else: #transient node
                types_sccs.append(0)
        types_sccs = np.array(types_sccs)
            
                
        
        dag = set()
        for ii,I in enumerate(Is_essential[i][:n_variables[i]]):
            for regulator in I:
                edge = (scc_dict[regulator],scc_dict[ii])
                if len(sccs[i])>40:
                    if edge[0]!=edge[1] and regulator<n_variables[i] and outdegrees[ii]>0:
                        dag.add(edge)
                else:
                    if edge[0]!=edge[1]:
                        dag.add(edge)                    
        nodes_in_dag = list(set.union(*[set(el) for el in dag]))
        
        labels = dict(zip(nodes_in_dag,[str(len_sccs[node]) if len_sccs[node]>1 else '' for node in nodes_in_dag]))
        g = nx.DiGraph(dag)
        top_to_bottom = [sorted(generation) for generation in nx.topological_generations(g)]
        pos = dict()
        for ii,generation in enumerate(nx.topological_generations(g)):
            len_generation = len(generation)
            if len_generation==1:
                pos.update( {generation[0] : (0,-ii)})
            else:
                delta_x = min(len_generation,7)/7
                for jj,node in enumerate(generation):
                    #pos.update( {node : (-0.5+jj/(len_generation-1),-ii)})
                    pos.update( {node : (-delta_x/2+delta_x*jj/(len_generation-1),-ii)})
        
        colors = np.array(['#555555','#ffaaaa','#555555','#ff7777'])[types_sccs[np.array(g.nodes)]+1]
        nx.draw_networkx(g,ax=ax[iii],labels=labels,pos=pos,node_color = colors,font_size=8,node_size=100*np.log(1.2*len_sccs[np.array(g.nodes)]))
        ax[iii].set_title(pmids_to_print[i].split(' ')[0]+r'$^*$'+' '.join(pmids_to_print[i].split(' ')[1:]) if len(sccs[i])>40 else pmids_to_print[i])
        ax[iii].axis('off')
    plt.savefig('all_sccs.pdf',bbox_inches = "tight")    


    
    
    
    
    
    
    
    
    
    
    
    









