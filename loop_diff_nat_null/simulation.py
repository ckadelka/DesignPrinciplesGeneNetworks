import load_database11 as db
import numpy as np
import canalizing_function_toolbox_v1_9 as can
import sys
import networkx as nx
import json

output_folder = 'results/'

max_loop = 6

def count_FBLs(Fs, Is, all_max_degrees, max_degree):
    ## Feedbackloops and triads
    all_loops = []
    all_types = []
    triads = []

    N = len(Fs)

    n_variables = [len(variables) for variables in variabless]
    for i in range(len(Fs)):
        F = Fs[i]
        I = Is[i]
        edges = []
        for j,regulators in enumerate(I):
            if j>=can.n_from_f(Fs[i]): #exclude constant self-loops
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
    return nr_real_loops.tolist()

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

def get_len_one(Is):
    M = can.adjacency_matrix(I)
    


folders = ['../update_rules_cell_collective/', '../update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 14
max_n = 40

Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded = db.load_database(folders,max_degree=max_degree,max_n=max_n)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)

Fs,Is,degrees = exclude_wlen0_funcs(Fs, Is, degrees)

all_max_degrees = [max(degree) for degree in degrees]


rndFs,rndIs = [],[]

filename = sys.argv[0]
SLURM_ID = int(sys.argv[1])

nsims_per_model = 1000
for i,F in enumerate(Fs):
    I = Is[i]
    n = [can.n_from_f(f) for f in F]
    n_essential = [can.nr_essential_variables(f) for f in F]
    k = [can.get_canalizing_depth(f) for f in F]
    
    rndFs.append([])
    rndIs.append([])
    for j in range(nsims_per_model):
        rndModel = can.random_BN(len(F), n = n_essential, k = k, STRONGLY_CONNECTED=False, NO_SELF_REGULATION=False,EXACT_DEPTH=True)
        rndFs[-1].append(rndModel[0])
        rndIs[-1].append(rndModel[1])
    print(str(100*(i/len(Fs)))+"%\t\t\t\t\t\t\t\t\t\t\t\t",end="\r")
    
FBLs = count_FBLs(Fs, Is, all_max_degrees, max_degree)
rndFBLs = []
for i,netFs in enumerate(rndFs):
    netIs = rndIs[i]
    rndFBLs.append(count_FBLs(netFs, netIs, all_max_degrees, max_degree))

f = open(output_folder+'SLURM_ID%i.txt' % SLURM_ID ,'w')
f.write('filename\t'+filename+'\n')
f.write('SLURM_ID\t'+str(SLURM_ID)+'\n')
vec = [json.dumps(FBLs), json.dumps(rndFBLs)]
name_vec = ['FBL counts','rndFBL counts']
for i in range(len(vec)):
    f.write(name_vec[i]+'\t'+vec[i]+'\n')
f.close()
