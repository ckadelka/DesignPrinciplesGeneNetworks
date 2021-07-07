import canalizing_function_toolbox_v1_9 as can
import numpy as np
import random as rand
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import sys
import time
import analysis_lib as analysis
import json

version = '02'
output_folder = 'results2/'


filename = sys.argv[0]
SLURM_ID = int(sys.argv[1])

print('Slurm ID: ' + str(SLURM_ID))
    
nsim = 1


#1) define parameters
N = 20
n = 2
k = n
nonModel = True

#Helper functions

def load_network():
    f = open('network.txt')
    s = f.read()
    f.close()
    return json.loads(s)

def record_loops(network):
    F = network[0]
    I = network[1]
    D = network[2]

    ret_length = []
    ret_type = []

    loops = analysis.get_loops(network)
    for loop in loops:
        l_length = len(loop)
        l_type = can.get_type_of_loop(loop,F,I)
        type_num = can.get_loop_type_number(l_type)
        ret_length.append(l_length)
        ret_type.append(type_num)
    return ret_length,ret_type

def record_attractors(F, I, N):
    ret = []
    attractors = can.num_of_attractors_v2(F, I, N)[0]
    for attractor in attractors:
        ret.append(len(attractor))
    return ret

def record_interactions(F):
    #have 10x3 array and record number of essential varaibles
    for f in F:
        t = can.is_monotonic2(f)
        interaction[0][can.nr_essential_variables(f)-1] += t[0]
        interaction[1][can.nr_essential_variables(f)-1] += t[1]
        interaction[2][can.nr_essential_variables(f)-1] += t[3]

#2) generate networks

net = load_network()

attractor_output = []
loop_output_length = []
loop_output_type = []
interaction = np.zeros(30).reshape(3,10)

def analyse_network(network):
    F = network[0]
    I = network[1]
    D = network[2]

    loop_len,loop_type = record_loops(network)
    loop_output_length.append(loop_len)
    loop_output_type.append(loop_type)

    attractor_output.append(record_attractors(F, I, len(F)))
    record_interactions(F)


TIME = time.time()
for i in range(nsim):
    if nonModel:
        for n in net:
            network = can.random_BN(n[0],n=n[1],k=n[2], STRONGLY_CONNECTED=False, NO_SELF_REGULATION=False, EXACT_DEPTH=True, MORE_ACTV=True)
            analyse_network(network)
    else:
        network = can.random_BN(N=N, n=n, k=k, indegree_distribution='constant')
        analyse_network(network)

print(interaction)
print(interaction/sum(interaction,0))

#2b) write results
f = open(output_folder+'results_v%s_nsim%i_SLURM_ID%i.txt' % (version,nsim,SLURM_ID) ,'w')
f.write('filename\t'+filename+'\n')
f.write('SLURM_ID\t'+str(SLURM_ID)+'\n')
f.write('nsim\t'+str(nsim)+'\n')
f.write('time in seconds\t'+str(int(time.time()-TIME))+'\n')
vec = [attractor_output, loop_output_type, loop_output_length, [list(a) for a in list(interaction)]]
name_vec = ['attractors','loop types','loop lengths', 'interactions']
for i in range(len(vec)):
    #f.write(name_vec[i]+'\t'+'\t'.join(list(map(str,vec[i])))+'\n')
    f.write(name_vec[i]+'\t'+json.dumps(vec[i])+'\n')
f.close()
