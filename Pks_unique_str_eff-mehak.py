#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 10:20:51 2022

@author: mkapoor
"""


import importlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import networkx as nx
import scipy.stats as stats
import canalizing_function_toolbox_v1_9 as can
import cana
import cana.boolean_node
import seaborn as sns
import numba as nb


n=4
bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
all_fs = np.array(list(itertools.product([0, 1], repeat=2**n)))
#all_fs = all_fs[:int(all_fs.shape[0]/2)]
strengths,Ps = [],[]
nr_symmetry_groups = []
depths = []
nr_essentials = []
res = []
eff_record =[]
for ii,f in enumerate(all_fs):
    dummy = can.get_canalizing_strength(f,bool_list=bool_list)
    strengths.append(dummy[0])
    nr_symmetry_groups.append(len(can.get_symmetry_groups(f,bool_list=bool_list)))
    Ps.append(dummy[1])
    depths.append(can.get_canalizing_depth_inputs_outputs_corefunction(f)[1])
    nr_essentials.append(can.nr_essential_variables(f))
    a = cana.boolean_node.BooleanNode(inputs=range(n),k=n,outputs=f)
    eff_rec = a.effective_connectivity()
    eff_record.append(eff_rec)
    layers=can.get_canalizing_depth_inputs_outputs_corefunction(f)
    #if ii%500==0:
        #print(ii)
    #if dummy[1][0]>0:
    #    res.append(sum(f))
    #    print(f,sum(f),nr_symmetry_groups[-1],dummy[1])
Ps = np.array(Ps)
strengths = np.array(strengths)
#print(np.mean(Ps[Ps[:,0]>0,:],0))
 
A = pd.DataFrame(np.c_[strengths,eff_record,Ps,all_fs.sum(1),nr_essentials,nr_symmetry_groups,depths,all_fs],columns=['strength']+['Effective connectivity']+list(['P_%i' % i for i in range(1,n)] +['Hamming weight','Nr essential variables','Nr symmetry groups','canalizing depth']+ [','.join(list(map(str,el))) for el in bool_list]))
A.to_csv('Get_canalizing_strength_of_all_Boolean_functions_with_n%i.csv' % n)
 
#extract all separate sequences of P_ks
#dummy = list(map(lambda x: 3*x[0]+5*x[1]+7*x[2],list(Ps)))
dummy_new = list(map(lambda x: 3*x[0]+5*x[1],list(Ps)))
dummy_pos = dict(zip(dummy_new,list(range(len(dummy_new)))))
already_included = [False]*len(dummy_new)
list_of_sequences = []
for i in range(2**(2**n)):
    if already_included[dummy_pos[dummy_new[i]]]==False:
        already_included[dummy_pos[dummy_new[i]]]=True
        list_of_sequences.append(Ps[i])
#A = pd.DataFrame(np.array(list_of_sequences))
#A.to_csv('Get_Pks_of_all_Boolean_functions_with_n%i_unique.csv' % n)
 
#extract one member of each equivalence class

dummy_new = list(map(lambda x: x[0]+3*x[1]+5*x[2]+7*x[3]+11*x[4],list(np.c_[Ps,depths,eff_record])))
dummy_pos = dict(zip(dummy_new,list(range(len(dummy_new)))))
already_included = [False]*len(dummy_new)
list_of_sequences = []
can_str=[]
can_eff=[]
for i in range(2**(2**n)):
    if already_included[dummy_pos[dummy_new[i]]]==False:
        already_included[dummy_pos[dummy_new[i]]]=True
        list_of_sequences.append([*Ps[i],strengths[i],eff_record[i],depths[i],*all_fs[i]])
        can_str.append(strengths[i])
        can_eff.append(eff_record[i])

A = pd.DataFrame(np.array(list_of_sequences),columns=list(['P_%i' % i for i in range(1,n)] +['strengths','effectiveness']+['canalizing depth']+ [','.join(list(map(str,el))) for el in bool_list]))
A.to_csv('Get_Pks_of_all_Boolean_functions_with_n%i_unique.csv' % n)

can_str=np.array(can_str)
can_eff=np.array(can_eff)


##Pareto Frontiers 
@nb.jit(nogil = True,cache = True,fastmath = True) 
def pareto(x,y,MAX_X=False,MAX_Y =False):
    indices = np.lexsort(((-1)**MAX_Y*y,(-1)**MAX_X*x))
    costs2 = [(-1)**MAX_Y*y[idx]for idx in indices]
    min_costs2 =costs2[0]
    pareto_indices =[indices[0]]
    for i,c in enumerate(costs2[1:]):
        if c<min_costs2:
            min_costs2=c
            pareto_indices.append(indices[i+1])
    return pareto_indices


##Call convex hull 
pareto_indices = pareto(can_str,can_eff,0,0)
print(pareto_indices)
f,ax = plt.subplots(figsize=(8,8))
ax.plot(can_str,can_eff,'ro')
ax.plot(can_str[pareto_indices],can_eff[pareto_indices],'bo')
ax.set_xlabel('Canalizing strength')
ax.set_ylabel('Effective Connectivity')
##Call Concave hull
pareto_indices2 = pareto(can_str,can_eff,1,1)
print(pareto_indices2)
f,ax = plt.subplots(figsize=(8,8))
ax.plot(can_str,can_eff,'ro')
ax.plot(can_str[pareto_indices2],can_eff[pareto_indices2],'go')
ax.set_xlabel('Canalizing strength')
ax.set_ylabel('Effective Connectivity')



