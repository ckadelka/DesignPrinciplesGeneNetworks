import canalizing_function_toolbox_v1_9 as can
import numpy as np
import random as rand
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import sys
import time
import analysis_lib as analysis

version = '01'
output_folder = 'results/'

try:
    filename = sys.argv[0]
    SLURM_ID = int(sys.argv[1])
except:
    filename = ''
    SLURM_ID = 1

if len(sys.argv)>2:
    nsim = int(sys.argv[2])
else:
    nsim = 100


#1) define parameters
N = 20
n = 6
k = n

#Helper functions

def record_loops(network):
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

def record_attractors(F, I):
    ret = []
    attractors = can.num_of_attractors_v2(F, I, N)[0]
    for attractor in attractors:
        ret.append(len(attractor))
    return ret

#2) generate networks

attractor_output = []
loop_output_length = []
loop_output_type = []

TIME = time.time()
for i in range(nsim):
    network = can.random_BN(N=N, n=n, k=k, indegree_distribution='constant')
    F = network[0]
    I = network[1]
    D = network[2]

    loop_len,loop_type = record_loops(network)
    loop_output_length.append(loop_len)
    loop_output_type.append(loop_type)

    attractor_output.append(record_attractors(F, I))



#2b) write results
f = open(output_folder+'results_v%s_nsim%i_SLURM_ID%i.txt' % (version,nsim,SLURM_ID) ,'w')
f.write('filename\t'+filename+'\n')
f.write('SLURM_ID\t'+str(SLURM_ID)+'\n')
f.write('nsim\t'+str(nsim)+'\n')
f.write('time in seconds\t'+str(int(time.time()-TIME))+'\n')
vec = [attractor_output, loop_output_type, loop_output_length]
name_vec = ['attractors','loop types','loop lengths']
for i in range(len(vec)):
    f.write(name_vec[i]+'\t'+'\t'.join(list(map(str,vec[i])))+'\n')
f.close()
