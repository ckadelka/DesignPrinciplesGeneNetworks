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
    
nsims = 5

folders = ['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 16
max_n= 1000
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded = db.load_database(folders,max_degree=max_degree,max_n=max_n)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)

print("loaded: ", N)
def analyze_networks(tFs, tIs, tdegrees):
	out = {
		"steady_states": [],
		"total_attractors": [],
		"nums_loop_type": [[] for i in range(5)]
	}
	for i in range(len(tFs)):
		F = tFs[i]
		I = tIs[i]
		D = tdegrees[i]
		network = [F, I, D]
		num_steady_state = 0
		num_attractors = 0
		num_loop_type = [0, 0, 0, 0, 0]

		length, loop_type = analysis.record_loops(network)
		for j in range(len(length)):
			if length[j] == 4 and loop_type[j] < 5:
				num_loop_type[loop_type[j]] += 1
		
		for j,num in enumerate(num_loop_type):
			out["nums_loop_type"][j].append(num)
		
		try:
			attractors = can.num_of_attractors_v2(F, I, len(F))
			for j in range(len(attractors[0])):
				if len(attractors[0][j]) == 1:
					num_steady_state += 1
			num_attractors += attractors[1]
			out["steady_states"].append(num_steady_state)
			out["total_attractors"].append(num_attractors)
		except:
			out["steady_states"].append(-1)
			out["total_attractors"].append(-1)
	
	return out


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

data = {
	"rewire": rewire(),
	"activation proportion": prop_pn(),
	"canalization": canalization()
}
f = open(output_folder+'results_v%s_nsim%i_SLURM_ID%i.txt' % (version,nsims,SLURM_ID) ,'w')
f.write(json.dumps(data))
f.close()
