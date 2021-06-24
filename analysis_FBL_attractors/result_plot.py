import os
import argparse
import ast
import numpy as np
import scipy.stats as st

import matplotlib
import matplotlib.pyplot as plt

#"python result_plot.py -d results -l 4 -t 0,1,2,3,4,5,6,7,8"
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

#x is 2d array [num loops][avg attr. sizes]
def get_x():
    x = [[] for i in range(22)]
    for k,i in enumerate(lengths):
        num_selected = i.count(selected_l)
        x[num_selected].append(np.mean(attrs[k]))
    return x

def loops_to_count():
    x = [[] for i in range(22)]
    for k,i in enumerate(lengths):
        num_selected = i.count(selected_l)
        x[num_selected].append(len(attrs[k]))
    return x

def scatterplot():
    f,ax = plt.subplots()
    ax.scatter(avg_attrs,loop,alpha=0.5,label='%i-loops' % (selected_l))
    ax.set_ylabel('%i-loops count' % (selected_l))
    ax.set_xlabel('avg attractor len')
    ax.set_yticks(range(np.max(loop)+1))
    plt.tight_layout()
    plt.savefig('total_avg_attr_%iloops_N%i.png' % (selected_l,len(avg_attrs)))

def boxplot(x):
    plt.boxplot(x)
    plt.ylabel('avg attractor len')
    plt.xlabel(str(selected_l) + ' loops count')
    plt.savefig('total_avg_attr_%iloops_N%i.png' % (selected_l,len(avg_attrs)), bbox_inches="tight")

def confidence_intervals(x):
    for i,num in enumerate(x):
        if len(num) < 3:
            continue
        try:
            if len(num) < 30:
                conf = st.t.interval(alpha=0.95, df=len(num)-1, loc=np.mean(num), scale=st.sem(num))
            else:
                conf = st.norm.interval(alpha=0.95, loc=np.mean(num), scale=st.sem(num))
            print(i, " ", conf, " ", len(num))
        except:
            print("failed on i for num ", num)



#violinplot
#boxplot()
#scatterplot()
x = loops_to_count()
confidence_intervals(x)