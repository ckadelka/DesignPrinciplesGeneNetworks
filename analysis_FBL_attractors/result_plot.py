import os
import argparse
import ast
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="Result Directory")
parser.add_argument("-l", help="Selected Length")
parser.add_argument("-t", help="Selected Types")
args = parser.parse_args()

attrs = []
types = []
lengths = []

for fname in os.listdir(args.d):
    if fname.endswith('.txt'):
        f = open(args.d+'/'+fname,'r')
        lines = f.readlines()
        for line in lines:
            if line.startswith('attractors'):	
                attrs += [ast.literal_eval(i) for i in line.split('\t')[1:]]
            if line.startswith("loop types"):	
                types += [ast.literal_eval(i) for i in line.split('\t')[1:]]
            if line.startswith("loop lengths"):	
                lengths += [ast.literal_eval(i) for i in line.split('\t')[1:]]

selected_l = int(args.l)
selected_t = [int(i) for i in args.t.split(',')]
loop = []
for i in range(len(lengths)):
    ls = lengths[i]
    cnt = 0
    for j in range(len(ls)):
        l = ls[j]
        t = types[i][j]
        if l == selected_l:
            if t in selected_t:
                cnt += 1
    loop.append(cnt)


avg_attrs = []
for attr in attrs:
    avg_attrs.append(np.mean(attr))

f,ax = plt.subplots()
ax.scatter(avg_attrs,loop,alpha=0.5,label='%i-loops' % (selected_l))
ax.set_ylabel('%i-loops count' % (selected_l))
ax.set_xlabel('avg attractor len')
ax.set_yticks(range(np.max(loop)+1))
plt.tight_layout()
plt.savefig('total_avg_attr_%iloops_N%i.pdf' % (selected_l,len(avg_attrs)))

