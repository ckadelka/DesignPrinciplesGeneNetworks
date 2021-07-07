import os
import argparse
import ast
import numpy as np
import scipy.stats as st
import json
import matplotlib
import matplotlib.pyplot as plt

#"python result_plot.py -d results -l 4 -t 0,1,2,3,4,5,6,7,8"
parser = argparse.ArgumentParser()
parser.add_argument("-d", help="Result Directory")
#parser.add_argument("-l", help="Selected Length")
#parser.add_argument("-t", help="Selected Types")
args = parser.parse_args()

attrs = []
types = []
lengths = []
interaction = np.zeros(30).reshape(3,10)

def calc_actv_inhb(interactions,name):
    n_max = 10
    color_unknown = 'black'#[0.5,0.5,0.5]
    color_neg = 'orange'#[0.7,0.7,1]
    color_pos = 'blue'#[1,0.7,0.7]
    colors = [color_pos,color_neg,color_unknown]
    f,ax = plt.subplots()
    x = np.arange(1,n_max+1)
    res_agg_prop = interactions/sum(interactions,0)
    for i,label in enumerate(['increasing','decreasing','not monotonic']):
        ax.bar(x,res_agg_prop[i,:],bottom=np.sum(res_agg_prop[:i,:],0),color=colors[i],alpha=0.3)
    ax.legend(['positive','negative','conditional'],loc=8,ncol=1,title='type of regulation')
    ax.set_ylim([0,1])
    ax.set_xticks(list(range(1,n_max+1)))
    ax.set_xlabel('number of essential inputs')
    ax.set_ylabel('proportion')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig(name,bbox_inches = "tight")

for fname in os.listdir(args.d):
    if fname.endswith('.txt'):
        f = open(args.d+'/'+fname,'r')
        lines = f.readlines()
        for line in lines:
            if line.startswith('attractors'):	
                attrs += json.loads(line[line.index('\t')+1:])
            if line.startswith("loop types"):	
                types += json.loads(line[line.index('\t')+1:])
            if line.startswith("loop lengths"):	
                lengths += json.loads(line[line.index('\t')+1:])
            if line.startswith("interactions"):	
                interact = json.loads(line[line.index('\t')+1:])
                interaction = np.add(interaction,interact)
print(len(lengths))
print(lengths)
print(len(types))
print(types)
calc_actv_inhb(interaction,'interactions.png')               
    
selected_l=4
loopAll = []
for i in range(len(lengths)):
    ls = lengths[i]
    cnt = 0
    for j in range(len(ls)):
        l = ls[j]
        t = types[i][j]
        if l == selected_l:
            loopAll.append(t)
loopAll2=[loopAll.count(i) for i in range(5)]
plt.clf()
plt.bar(range(5),loopAll2)
plt.ylabel('4 loop count')
plt.xlabel('4 loop type')
plt.savefig('4loopTypes.png')

