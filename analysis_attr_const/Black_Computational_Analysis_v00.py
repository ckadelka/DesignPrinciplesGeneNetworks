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

folders = ['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 16
max_n=12
[Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded] = db.load_database(folders,max_degree=max_degree,max_n=max_n)

#Rewire I
rewired_steady_states = []
rewired_total_attractors = []
rewired_num_loop_type_0 = []
rewired_num_loop_type_4 = []
rewired_num_loop_type_1 = []
rewired_num_loop_type_3 = []
for i in range(len(Fs)):
	F = Fs[i]
	I = can.rewire_I(Is[i])
	D = degrees[i]
	network = [F, I, D]
	num_steady_state = 0
	num_attractors = 0
	num_loop_type_0 = 0
	num_loop_type_4 = 0
	num_loop_type_1 = 0
	num_loop_type_3 = 0

	length, loop_type = analysis.record_loops(network)
	for j in range(len(length)):
		if length[j] == 4 and loop_type[j] == 0:
			num_loop_type_0 += 1
		if length[j] == 4 and loop_type[j] == 4:
			num_loop_type_4 += 1
		if length[j] == 4 and loop_type[j] == 1:
			num_loop_type_1 += 1
		if length[j] == 4 and loop_type[j] == 3:
			num_loop_type_3 += 1
	rewired_num_loop_type_0.append(num_loop_type_0)
	rewired_num_loop_type_4.append(num_loop_type_4)
	rewired_num_loop_type_1.append(num_loop_type_1)
	rewired_num_loop_type_3.append(num_loop_type_3)

	attractors = can.num_of_attractors_v2(F, I, N)
	for j in range(len(attractors[0])):
		if len(attractors[0][j]) == 1:
			num_steady_state += 1
	num_attractors += attractors[1]
	rewired_steady_states.append(num_steady_state)
	rewired_total_attractors.append(num_attractors)


#Generate networks preserving size, indegree distribution, k, and proportion +/- interactions
interactions_preserved_steady_states = []
interactions_preserved_total_attractors = []
interactions_preserved_num_loop_type_0 = []
interactions_preserved_num_loop_type_4 = []
interactions_preserved_num_loop_type_1 = []
interactions_preserved_num_loop_type_3 = []
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
	initial_proportion = increasing / decreasing
	network = can.random_BN(N = len(Fs[i]), n = degrees[i], k = [can.get_canalizing_depth(f) for f in Fs[i]], indegree_distribution = 'poisson', activation_proportion = initial_proportion)
	F = network[0]
	I = network[1]
	D = network[2]
	num_steady_state = 0
	num_attractors = 0
	num_loop_type_0 = 0
	num_loop_type_4 = 0
	num_loop_type_1 = 0
	num_loop_type_3 = 0

	length, loop_type = analysis.record_loops(network)
	for j in range(len(length)):
		if length[j] == 4 and loop_type[j] == 0:
			num_loop_type_0 += 1
		if length[j] == 4 and loop_type[j] == 4:
			num_loop_type_4 += 1
		if length[j] == 4 and loop_type[j] == 1:
			num_loop_type_1 += 1
		if length[j] == 4 and loop_type[j] == 3:
			num_loop_type_3 += 1
	interactions_preserved_num_loop_type_0.append(num_loop_type_0)
	interactions_preserved_num_loop_type_4.append(num_loop_type_4)
	interactions_preserved_num_loop_type_1.append(num_loop_type_1)
	interactions_preserved_num_loop_type_3.append(num_loop_type_3)

	attractors = can.num_of_attractors_v2(F, I, N)
	for j in range(len(attractors[0])):
		if len(attractors[0][j]) == 1:
			num_steady_state += 1
	num_attractors += attractors[1]
	interactions_preserved_steady_states.append(num_steady_state)
	interactions_presrved_total_attractors.append(num_attractors)

#Generate networks preserving size, indegree distribution, and k
k_preserved_steady_states = []
k_preserved_total_attractors = []
k_preserved_num_loop_type_0 = []
k_preserved_num_loop_type_4 = []
k_preserved_num_loop_type_1 = []
k_preserved_num_loop_type_3 = []
for i in range(len(Fs)):
	network = can.random_BN(N = len(Fs[i]), n = degrees[i], k = [can.get_canalizing_depth(f) for f in Fs[i]], indegree_distribution = 'poisson')
	F = network[0]
	I = network[1]
	D = network[2]
	num_steady_state = 0
	num_attractors = 0
	num_loop_type_0 = 0
	num_loop_type_4 = 0
	num_loop_type_1 = 0
	num_loop_type_3 = 0

	length, loop_type = analysis.record_loops(network)
	for j in range(len(length)):
		if length[j] == 4 and loop_type[j] == 0:
			num_loop_type_0 += 1
		if length[j] == 4 and loop_type[j] == 4:
			num_loop_type_4 += 1
		if length[j] == 4 and loop_type[j] == 1:
			num_loop_type_1 += 1
		if length[j] == 4 and loop_type[j] == 3:
			num_loop_type_3 += 1
	k_preserved_num_loop_type_0.append(num_loop_type_0)
	k_preserved_num_loop_type_4.append(num_loop_type_4)
	k_preserved_num_loop_type_1.append(num_loop_type_1)
	k_preserved_num_loop_type_3.append(num_loop_type_3)

	attractors = can.num_of_attractors_v2(F, I, N)
	for j in range(len(attractors[0])):
		if len(attractors[0][j]) == 1:
			num_steady_state += 1
	num_attractors += attractors[1]
	k_preserved_steady_states.append(num_steady_state)
	k_preserved_total_attractors.append(num_attractors)
