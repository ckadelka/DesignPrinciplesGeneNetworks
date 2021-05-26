#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:01:32 2019

@author: ckadelka
"""

##Imports

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import canalizing_function_toolbox_v1_9 as can
import load_database09 as db
import itertools

plt.rcParams.update({'font.size': 16})

## load the database, choose low max_n for quick results and to only look at small models
folders = ['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 16
max_n=1000
[Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded] = db.load_database(folders,max_degree=max_degree,max_n=max_n)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)


## calculate the "degree" of canalization for each function
depths,layer_structures,textfiles,all_ncfs = [],[],[],[]
ks_per_n = np.zeros((max_degree+15,max_degree+15))
for i in range(N):
    constants = constantss[i]
    variables = variabless[i]
    F = Fs[i]
    depths.append([])
    layer_structures.append([])
    for f in F:
        if len(f)==0: #happens if actual degree_f > max_degree
            depth = np.nan
            layer_structure = []
        else:
            (n_f,depth,can_inputs,can_outputs,corefunction) = can.get_canalizing_depth_inputs_outputs_corefunction(f)
            layer_structure = can.get_layer_structure_given_outputs_corefunction(can_outputs,corefunction,n_f)
        depths[-1].append(depth)
        layer_structures[-1].append(layer_structure)
    all_ncfs.append(sum(depths[-1]) == sum(degrees[-1]))
    for k,n in zip(depths[-1],degrees_essential[i]):
        if not np.isnan(k):
            ks_per_n[n,k] += 1

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


upto = 10
nsim = 1000
ks_per_n_random = np.zeros((upto+1,upto+1))
ss_per_n_random = np.zeros((upto+1,upto+1))
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

res = []
for ii,F in enumerate(Fs):
    for f in F:
        (n,k,can_inputs,can_outputs,corefunction) = can.get_canalizing_depth_inputs_outputs_corefunction(f)
        if k>-1:
            kis = can.get_layer_structure_given_outputs_corefunction(can_outputs,corefunction,n)
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
        res_rnd.append(can.get_layer_structure_given_outputs_corefunction(can.get_canalizing_depth_inputs_outputs_corefunction(can.random_k_canalizing(nk,nk,True))[3],[0],nk))
    a = pd.Series(res_rnd)
    print(a.value_counts())
    
    a = pd.Series(list(map(len,res_rnd)))
    print(a.value_counts())
    dummy = a.value_counts()
    for key in dummy.keys():
        exp_counts[nk-3,key-1]  = dummy[key]
        
pd.DataFrame(obs_counts,index=['n = '+el for el in list(map(str,range(3,upto+1)))],columns=list(map(str,range(1,upto+1)))).to_excel('observed_number_of_layers_n%i.xlsx' % len(Fs))
pd.DataFrame(exp_counts,index=['n = '+el for el in list(map(str,range(3,upto+1)))],columns=list(map(str,range(1,upto+1)))).to_excel('expected_number_of_layers_nsim%i.xlsx' % (nsim))



nk = 3
a = pd.Series([list(el[:-1]) for el in res[np.bitwise_and(res[:,1]==nk,res[:,2]==nk)][:,4]])
ex = a.value_counts()
pd.DataFrame(np.c_[list(ex.index),ex.values],columns = ['b'+str(i+1) for i in range(nk-1)]+['count']).to_excel('can_output_distribution_for_ncfs_with_n%i.xlsx' % nk)

nk = 3
a = pd.Series([list(el[:-1]) for el in res[np.bitwise_and(res[:,1]==nk,res[:,2]==nk)][:,3]])
ex = a.value_counts()
pd.DataFrame(np.c_[list(ex.index),ex.values],columns = ['a'+str(i+1) for i in range(nk-1)]+['count']).to_excel('can_input_distribution_for_ncfs_with_n%i.xlsx' % nk)

f,ax = plt.subplots()
for nk in range(1,11):
    ax.bar([nk],[np.mean(res[np.bitwise_and(res[:,1]==nk,res[:,2]==nk)][:,7])],color='k')
ax.set_xlabel('number of variables')
ax.set_ylabel('average proportion of 1s')
plt.gcf().subplots_adjust(bottom=0.2)
plt.savefig('avg1s_vs_n.pdf')
    
    
    





#calculate symmetry groups
nr_symmetry_groups =[]
upto=11
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

ss_per_n_random = np.zeros((upto+1,upto+1))
nsim = 1000
for n in range(1,upto+1):
    bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
    for _ in range(nsim):
        F = can.random_non_degenerated_function(n)
        s = len(can.get_symmetry_groups(F,bool_list))
        ss_per_n_random[n,s] += 1
pd.DataFrame(ss_per_n_random,index=['n=%i' % n for n in range(len(ss_per_n_random))],columns=['s=%i' % n for n in range(len(ss_per_n_random))]).to_excel('nr_symmetry_groups_per_n_random_nsim%i.xlsx' % (nsim))

ss_per_n_random_canalizing = np.zeros((upto+1,upto+1))
nsim = 1000
for n in range(1,upto+1):
    bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
    for _ in range(nsim):
        F = can.random_k_canalizing(n,1,False,bool_list)
        s = len(can.get_symmetry_groups(F,bool_list))
        ss_per_n_random_canalizing[n,s] += 1
pd.DataFrame(ss_per_n_random_canalizing,index=['n=%i' % n for n in range(len(ss_per_n_random_canalizing))],columns=['s=%i' % n for n in range(len(ss_per_n_random_canalizing))]).to_excel('nr_symmetry_groups_per_n_random_canalizing_nsim%i.xlsx' % (nsim))

ss_per_n_random_NCF = np.zeros((upto+1,upto+1))
nsim = 1000
for n in range(1,upto+1):
    bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
    for _ in range(nsim):
        F = can.random_k_canalizing(n,n,False,bool_list)
        s = len(can.get_symmetry_groups(F,bool_list))
        ss_per_n_random_NCF[n,s] += 1
pd.DataFrame(ss_per_n_random_NCF,index=['n=%i' % n for n in range(len(ss_per_n_random_NCF))],columns=['s=%i' % n for n in range(len(ss_per_n_random_NCF))]).to_excel('nr_symmetry_groups_per_n_random_NCF_nsim%i.xlsx' % (nsim))

ss_per_n_imputed = np.zeros((upto+1,upto+1))
nsim = 1000
for n in range(1,upto+1):
    bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
    for k in np.random.choice(list(range(n+1)),nsim,replace=True,p=ks_per_n[n][:(n+1)]/sum(ks_per_n[n][:(n+1)])):
        F = can.random_k_canalizing(n,k,True,bool_list)
        s = len(can.get_symmetry_groups(F,bool_list))
        ss_per_n_imputed[n,s] += 1
pd.DataFrame(ss_per_n_imputed,index=['n=%i' % n for n in range(len(ss_per_n_imputed))],columns=['s=%i' % n for n in range(len(ss_per_n_imputed))]).to_excel('nr_symmetry_groups_per_n_random_imputed_nsim%i.xlsx' % (nsim))






## Canalizing strength - needs work!  
can_strengths = []
for F in Fs:
    can_strengths.append([])
    for f in F:
        len_f = len(f)
        if len_f>2 and len_f<33:
            can_strengths[-1].append(can.get_canalizing_strength(f))

N=100
can_strengths = []
for n in range(2,6):
    can_strengths.append([])
    for _ in range(N):
        can_strengths[-1].append( can.get_canalizing_strength(can.random_k_canalizing(n,1)) )

can_strengths_single = []
for i in range(len(range(2,11))):
    can_strengths_single.append([el[0] for el in can_strengths[i]])


#canalizing strength of non-canalizing functions (n>=3)
a,b = 3,6
nsim = 100
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


import matplotlib.patches as mpatches
f,ax = plt.subplots()

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
ax.legend(*zip(*labels),loc='best')
plt.gcf().subplots_adjust(bottom=0.17,left=0.2)
ax.set_title('')
plt.savefig('canalizing_strength_nmin%i_nmax%i_k%i_N%i_nsim%i.pdf' % (a,b,k,N,nsim))







