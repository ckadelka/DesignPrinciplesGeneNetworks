import canalizing_function_toolbox_v1_9 as can
import analysis_lib as analysis
import load_database11 as db
import numpy as np
import random as rand
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib.ticker import PercentFormatter
from matplotlib import cm
import sys
import time
import networkx as nx
import json

version = '01'
output_folder = 'results/'


filename = sys.argv[0]
SLURM_ID = int(sys.argv[1])
TIME = time.time()
    
nsims = 1

folders = ['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 12
max_n= 30
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded = db.load_database(folders,max_degree=max_degree,max_n=max_n)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
# Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)

print("loaded: ", N)
def analyze_networks(tFs, tIs, tdegrees):
	out = {
		"steady_states": [],
		"total_attractors": [],
		"avg_length_attractors": [],
		"entropy": [],
		"nums_loop_type": [[] for i in range(5)],
		"all_ffls": []
	}
	for i in range(len(tFs)):
		F = tFs[i]
		I = tIs[i]
		D = tdegrees[i]
		network = [F, I, D]
		num_steady_state = 0
		num_attractors = 0
		avg_len_attractors = 0
		total_attractor_length = 0
		num_loop_type = [0, 0, 0, 0, 0]
		

		length, loop_type = analysis.record_loops(network)
		for j in range(len(length)):
			if length[j] == 4 and loop_type[j] < 5:
				num_loop_type[loop_type[j]] += 1
		
		for j,num in enumerate(num_loop_type):
			out["nums_loop_type"][j].append(num)
		
		for i in range(len(Fs)):
			F = Fs[i]
			I = Is[i]
			A = can.adjacency_matrix(I,constantss[i])
			(ffls,types) = can.get_ffls(A,F,I)
			out["all_ffls"].append(list(map(can.get_ffl_type_number,types)))
		
		
		try:
			attractors = can.num_of_attractors_v2(F, I, len(F))
			for j in range(len(attractors[0])):
				if len(attractors[0][j]) == 1:
					num_steady_state += 1
			num_attractors += attractors[1]
			out["steady_states"].append(num_steady_state)
			out["total_attractors"].append(num_attractors)

			for j in range(attractors[1]):
				total_attractor_length += len(attractors[0][j])
			out["avg_length_attractors"].append(total_attractor_length / attractors[1])

			out["entropy"].append(can.entropy(attractors[2]))
		except:
			out["steady_states"].append(-1)
			out["total_attractors"].append(-1)
			out["avg_length_attractors"].append(-1)
			out["entropy"].append(-1)
		
	return out

def real_networks():
	return analyze_networks(Fs, Is, degrees)

#rewire I
def rewire():
	nFs = []
	nIs = []
	ndegrees = []
	for i in range(nsims):
		for i,I in enumerate(Is):
			nIs.append(can.rewire_I(I,preserve_self_regulation=True,NO_SELF_REGULATION=True))
			nFs.append(Fs[i])
			ndegrees.append(degrees[i])
	return analyze_networks(Fs, nIs, degrees)
		

#Generate networks preserving size, indegree distribution, k, and proportion +/- interactions
def prop_pn():
	tFs = []
	tIs = []
	tdegrees = []
	for j in range(nsims):
		for i in range(len(Fs)):
			types = [can.variable_types(f) for f in Fs[i]]
			increasing = 0
			decreasing = 0
			for j in range(len(types)):
				for k in range(len(types[j])):
					if types[j][k] == 'decreasing':
						decreasing += 1
					elif types[j][k] == 'increasing':
						increasing += 1
			initial_proportion = increasing / (increasing + decreasing)
			network = can.random_BN(N = len(Fs[i]), n = degrees[i], k = [can.get_canalizing_depth(f) for f in Fs[i]], indegree_distribution = 'poisson', activation_proportion = initial_proportion)
			F = network[0]
			I = network[1]
			D = network[2]
			tFs.append(F)
			tIs.append(I)
			tdegrees.append(D)
	return analyze_networks(tFs, tIs, tdegrees)

#Generate networks preserving size, indegree distribution, and k
def canalization():
	tFs = []
	tIs = []
	tdegrees = []
	for j in range(nsims):
		for i in range(len(Fs)):
			network = can.random_BN(N = len(Fs[i]), n = degrees[i], k = [can.get_canalizing_depth(f) for f in Fs[i]], indegree_distribution = 'poisson')
			tFs.append(network[0])
			tIs.append(network[1])
			tdegrees.append(network[2])

	return analyze_networks(tFs, tIs, tdegrees)

data = [rewire(), prop_pn(), canalization()]
data_names = ["rewire", "activation proportion", "canalization"]

for i,d in enumerate(data):
	data_name = data_names[i]
	f = open(output_folder+'results_operation-'+data_name+'_nsim%i_SLURM_ID%i.txt' % (nsims,SLURM_ID) ,'w')
	f.write('filename\t'+filename+'\n')
	f.write('SLURM_ID\t'+str(SLURM_ID)+'\n')
	f.write('nsim\t'+str(nsims)+'\n')
	f.write('time in seconds\t'+str(int(time.time()-TIME))+'\n')

	name_vec = [key for key in d.keys()]
	vec = [d[key] for key in name_vec]

	for i in range(len(vec)):
		f.write(name_vec[i]+'\t'+'\t'.join(list(map(str,vec[i])))+'\n')
	f.close()
