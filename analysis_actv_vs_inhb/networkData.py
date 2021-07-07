import load_database11 as db
import numpy as np
import canalizing_function_toolbox_v1_9 as can
import sys
import networkx as nx
import json

def exclude_wlen0_funcs(Fs,Is,degrees):
    oFs = []
    oIs = []
    oDegrees = []
    for i,F in enumerate(Fs):
        I = Is[i]
        n_essential = [can.nr_essential_variables(f) for f in F]
        if min(n_essential) > 0:
            oFs.append(F)
            oIs.append(I)
            oDegrees.append(degrees[i])
    return oFs,oIs,oDegrees


folders = ['../update_rules_cell_collective/', '../update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 14
max_n = 40

Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded = db.load_database(folders,max_degree=max_degree,max_n=max_n)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)

Fs,Is,degrees = exclude_wlen0_funcs(Fs, Is, degrees)

all_max_degrees = [max(degree) for degree in degrees]
print('all_max_degrees '  + str(all_max_degrees))
maxD=[]
for F in Fs:
    for f in F:
        maxD.append(can.n_from_f(f))
print('Max degree: ' + str(max(maxD)))


res = []
for i,F in enumerate(Fs):
    I = Is[i]
    n = [can.n_from_f(f) for f in F]
    n_essential = [can.nr_essential_variables(f) for f in F]
    k = [can.get_canalizing_depth(f) for f in F]
    res.append([len(F), n_essential, k])

f = open('network.txt','w')
f.write(json.dumps(res))
f.close()
