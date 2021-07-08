import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import canalizing_function_toolbox_v1_9 as can
import networkx as nx
import load_database11 as db

## load the database, choose low max_n for quick results and to only look at small models
#folders = ['mesih/']
folders = ['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 16
max_n=1000
[Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded] = db.load_database(folders,max_degree=max_degree,max_n=max_n,DONT_LOAD_IFMAX_DEGREE=True)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)

## Feedbackloops and triads
max_loop=4
all_loops = []
all_types = []
triads = []
n_variables = [len(variables) for variables in variabless]
all_max_degrees = [max(degree) for degree in degrees]
for i in range(len(Fs)):
    if all_max_degrees[i]>max_degree:
        all_types.append([])
        all_loops.append([])
        continue
    F = Fs[i]
    I = Is[i]
    degree = degrees[i]
    edges = []
    for j,regulators in enumerate(I):
        if j>=n_variables[i]: #exclude constant self-loops
            break
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

attr = []
for i in range(len(Fs)):
    F = Fs[i]
    I = Is[i]
    attr.append(len(can.num_of_attractors_v2(F,I,len(F))[0]))

for k in range(max_loop): 
    f,ax = plt.subplots()
    ax.scatter(attr,nr_loops[k],alpha=0.5,label='%i-loops' % (k+1))
    ax.set_ylabel('avg attr size in %i-loops' % (k+1))
    plt.tight_layout()
    plt.savefig('total_avg_attr_%iloops_N%i_v2.pdf' % (k+1,N))

import scipy.stats as stats
vec = [attr,nr_loops[3]]
len_vec = len(vec)
spearman = np.ones((len_vec,len_vec))
for i in range(len_vec):
    for j in range(i+1,len_vec):
        spearman[i,j] = stats.spearmanr(vec[i],vec[j])[0]
        spearman[j,i] = spearman[i,j]   

print(spearman)
