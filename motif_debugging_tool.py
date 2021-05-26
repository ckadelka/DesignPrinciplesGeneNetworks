# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:14:25 2020

@author: tbutrie
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import canalizing_function_toolbox_v1_9 as can
import networkx as nx
import load_database09 as db

file_name = input('Enter filename here ') + '.txt'
#Potential folder names 
#update_rules_cell_collective/ update_rules_models_in_literature_we_randomly_come_across/
folder_name = input('Enter folder name ')
F, I, degree, variables, constants = db.text_to_BN(folder_name, file_name)
G=can.generate_networkx_graph_from_edges(I,len(variables))
nx.draw(G)

max_loop=4
all_loops = []
all_types = []
triads = []
edges = []
for j,regulators in enumerate(I):
    for ii in regulators:
        edges.append((ii,j))
G=nx.DiGraph(edges)
triads.append(nx.triadic_census(G))
loops = list(can.simple_cycles(G,max_loop))
all_types.append([])
for loop in loops:
    all_types[-1].append(can.get_type_of_loop(loop,F,I))
all_loops.append(loops)
    
all_triads_keys = triads[0].keys()
all_triads_counts = [sum([el[key] for el in triads]) for key in all_triads_keys]

## Clustering of two FFLs or a FFL and a 2-loop or a FFL and a 3-loop
all_ffls = []
type_ffl_ffls = []
type_ffl_2l = []
type_ffl_3l = []

DEBUG = True

A = can.adjacency_matrix(I,constants)
(ffls,types) = can.get_ffls(A,F,I)

type_ffl_ffl = []
all_ffls.append(list(map(can.get_ffl_type_number,types)))
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

#FFLs clustered with 3 loops
max_loop_length = 3

loops = list(can.simple_cycles(G,max_len = max_loop_length))
two_loops = [el for el in loops if len(el)==2]
three_loops = [el for el in loops if len(el)==3]

for i in range(n_ffls):
    for j in range(len(two_loops)):            
        intersect = set.intersection(set(ffls[i]),set(two_loops[j]))
        if len( intersect )>0:
            types = dict(zip(ffls[i],['r','i','t']))
            for el in list(intersect):
                types.update({el:'m'})
            sorted_types = list(types.values())
            sorted_types.sort()
            string_type = ''.join(sorted_types) 
            print(string_type)
            if string_type == 'mmr':
                type_ffl_2l.append(13)
            elif string_type == 'mmt':
                type_ffl_2l.append(14)
            elif string_type == 'imm':
                type_ffl_2l.append(18)
            elif string_type == 'imt':
                type_ffl_2l.append(19)
            elif string_type == 'mrt':
                type_ffl_2l.append(20)
            elif string_type == 'imr':
                type_ffl_2l.append(21)
#Comparing Triads to Motif numbers
# 13 = 120D
# 14 = 120U
# 15/16 = 210
# 17 = 300
# 18 = 120C
#pd.value_counts()
#list(zip(all_triads_keys,all_triads_counts))

#A = np.c_[list(all_triads_keys),all_triads_counts]
#pd.DataFrame(A,columns=['key','count']).to_excel('fsdf.xlsx')