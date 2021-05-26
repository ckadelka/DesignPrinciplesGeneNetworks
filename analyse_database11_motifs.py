#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:01:32 2019

@author: ckadelka
"""

##Imports

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import canalizing_function_toolbox_v1_9 as can
import networkx as nx
import load_database11 as db

plt.rcParams.update({'font.size': 16})

## load the database, choose low max_n for quick results and to only look at small models
folders = ['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 16
max_n=1000
[Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded] = db.load_database(folders,max_degree=max_degree,max_n=max_n)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)

n_variables = [len(variables) for variables in variabless]
all_avg_degrees = [np.mean(degree) for degree in degrees]
all_max_degrees = [max(degree) for degree in degrees]

## find all FFLs
all_ffls = []
type_ffl = []
for i in range(len(Fs)):
    if all_max_degrees[i]>max_degree:
        all_ffls.append([])
        continue
    F = Fs[i]
    I = Is[i]
    A = can.adjacency_matrix(I,constantss[i])
    (ffls,types) = can.get_ffls(A,F,I)
    all_ffls.append(list(map(can.get_ffl_type_number,types)))

all_ffls_flat = []
for el in all_ffls:
    all_ffls_flat.extend(el)
    
LEGEND = list(itertools.product(['decreasing', 'increasing'], repeat=3))
LEGEND_NR = list(map(can.get_ffl_type_number,LEGEND))
LEGEND_COH = np.array(list(map(can.is_ffl_coherent,LEGEND)))

nr_ffls = np.array([len(el) for el in all_ffls])
nr_coh_ffls = np.zeros(N,dtype=int)
nr_incoh_ffls = np.zeros(N,dtype=int)
nr_unkown_ffls = np.zeros(N,dtype=int)
nr_notreal_ffls = np.zeros(N,dtype=int)
nr_specific_ffl = np.zeros((8,N),dtype=int)
for ii,ffls in enumerate(all_ffls):
    for el in ffls:
        if el==-1:
            nr_unkown_ffls[ii]+=1
        elif el==-2:
            nr_notreal_ffls[ii]+=1
        elif LEGEND_COH[el] == True:
            nr_coh_ffls[ii] += 1
            nr_specific_ffl[el,ii] += 1
        else:
            nr_incoh_ffls[ii] += 1
            nr_specific_ffl[el,ii] += 1
nr_real_ffls = nr_ffls-nr_notreal_ffls



## plotting FFLs
#NR FFL vs Network Size
f,ax = plt.subplots()
ax.scatter(n_variables,nr_real_ffls,alpha=0.5)

#NR FFL vs Average degree
f,ax = plt.subplots()
ax.scatter(all_avg_degrees,nr_real_ffls,alpha=0.5)

#count coherent vs incoherent vs unknown stratified by model
color_unknown = [0.5,0.5,0.5]
color_incoh = [0.7,0.7,1]
color_coh = [1,0.7,0.7]

DONT_SHOW_ZERO_FFLS_NETWORKS = True
LOG=False
sorted_sizes = np.array(sorted(zip(nr_real_ffls,nr_coh_ffls,nr_incoh_ffls,nr_unkown_ffls),key=lambda x: (x[0],x[1]))).T
if DONT_SHOW_ZERO_FFLS_NETWORKS:
    index = list(sorted_sizes[0]>0).index(True)
    sorted_sizes = sorted_sizes[:,index:]
f,ax=plt.subplots(figsize=(10,4))
ind = np.arange(sorted_sizes.shape[1])
ax.bar(ind,sorted_sizes[1],width=0.7,color=color_coh,label='coherent FFLs',log=LOG)
ax.bar(ind,sorted_sizes[2],width=0.7,bottom=sorted_sizes[1],color=color_incoh,label='incoherent FFLs',log=LOG)
ax.bar(ind,sorted_sizes[3],width=0.7,bottom=sorted_sizes[1]+sorted_sizes[2],color=color_unknown,label='unknown type FFL',log=LOG)
ax.set_xlim([-1,sorted_sizes.shape[1]])
ax.set_ylabel('Number of FFLs')
ax.set_xticks([])
ax.set_xlabel('Gene regulatory networks')
ax.legend(loc=2)
plt.savefig('number_of_coherent_ffls_N%i.pdf' % N)

#proportion coherent vs incoherent vs unknown stratified by model
DONT_SHOW_ZERO_FFLS_NETWORKS = True
LOG=False
nr_ffls-nr_notreal_ffls
sorted_sizes = np.array(sorted(zip(nr_real_ffls,nr_coh_ffls/nr_real_ffls,nr_incoh_ffls/nr_real_ffls,nr_unkown_ffls/nr_real_ffls),key=lambda x: (x[0],x[1]))).T
if DONT_SHOW_ZERO_FFLS_NETWORKS:
    index = list(sorted_sizes[0]>0).index(True)
    sorted_sizes = sorted_sizes[:,index:]
f,ax=plt.subplots(figsize=(10,4))
ind = np.arange(sorted_sizes.shape[1])
ax.bar(ind,sorted_sizes[1],width=0.7,color=color_coh,label='coherent FFLs',log=LOG)
ax.bar(ind,sorted_sizes[2],width=0.7,bottom=sorted_sizes[1],color=color_incoh,label='incoherent FFLs',log=LOG)
ax.bar(ind,sorted_sizes[3],width=0.7,bottom=sorted_sizes[1]+sorted_sizes[2],color=color_unknown,label='unknown type FFL',log=LOG)
ax.set_xlim([-1,sorted_sizes.shape[1]])
ax.set_ylabel('Proportion of FFLs')
ax.set_xticks([])
ax.set_xlabel('Gene regulatory networks')
ax2 = ax.twinx()
ax2.semilogy(ind,sorted_sizes[0],'x-',lw=2,color=[0.2,0.9,0.2])
ax2.set_ylabel('Number of FFLs',color=[0.2,0.9,0.2])
ax2.tick_params(axis='y', labelcolor=[0.2,0.9,0.2])
ax.legend(loc=8)
plt.savefig('proportion_of_coherent_ffls_N%i.pdf' % N)

#total count, coherent vs incoherent vs unknown
f,ax = plt.subplots()
ax.bar([0],[sum(nr_coh_ffls)],color=color_coh)
ax.bar([1],[sum(nr_incoh_ffls)],color=color_incoh)
ax.bar([2],[sum(nr_unkown_ffls)],color=color_unknown)
ax.set_xticks(range(3))
ax.set_xticklabels(['coherent','incoherent','unknown'],rotation=90)
ax.set_ylabel('total count')
plt.savefig('total_count_of_ffls_N%i.pdf' % N)

#total count, specific FFL type
order = np.append(np.arange(8)[LEGEND_COH],np.arange(8)[~LEGEND_COH])
cmap = matplotlib.cm.Paired

f,ax = plt.subplots()
for i in range(8):
    ax.bar([i],[sum(nr_specific_ffl[order[i],:])],color=cmap(i))
ax.bar([8],[sum(nr_unkown_ffls)],color=color_unknown)
ax.set_xticks(range(9))
ax.set_xticklabels(['type '+str(order[i]+1) for i in range(8)]+['unknown'],rotation=90)
ax.set_ylabel('total count')
plt.savefig('total_count_of_specific_ffls_N%i.pdf' % N)

#total count, specific FFL type - nice
import matplotlib
def arrow_new(self, x, y, dx, dy, **kwargs):
    kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
    kwargs.setdefault('fc', 'black')
    x = self.convert_xunits(x)
    y = self.convert_yunits(y)
    dx = self.convert_xunits(dx)
    dy = self.convert_yunits(dy)
    posA = x, y
    posB = x+dx, y+dy
    a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
    self.add_artist(a)
    return a

f,ax = plt.subplots()
width=0.8
for i in range(8):
    ax.bar([i],[sum(nr_specific_ffl[order[i],:])],color=cmap(i),width=width)
ax.bar([8],[sum(nr_unkown_ffls)],color=color_unknown)
ax.set_xticks(range(9))
ax.set_xticklabels(['' for i in range(8)]+['unknown'],rotation=90)
ax.xaxis.set_ticks_position('none') 
ax.set_ylabel('total count')
ax2 = ax.twinx()
yoffset = 0.13
epsilon = 0.25
total_length_y = 0.2
activation_head_length = 7
inhibition_head_length = 0.1
ax2.plot([0-width/2,3+width/2],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
ax2.plot([4-width/2,7+width/2],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
ax2.text(1.5,1.45*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
ax2.text(5.5,1.45*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
ax2.set_ylim([0,1])
for j in range(8):
    ax2.text(j,-.07,str(j+1),va='center',ha='center',color=cmap(j),clip_on=False)
for i in range(8):
    color=cmap(i)
    direct,indirect1,indirect2 = LEGEND[order[i]]
    head_width = 7
    head_length = activation_head_length if direct == 'increasing' else inhibition_head_length
    arrow_new(ax2,i+epsilon,-yoffset,0,-total_length_y,clip_on=False,fc=color,ec=color,color=color,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    head_length = activation_head_length if indirect1 == 'increasing' else inhibition_head_length
    arrow_new(ax2,i+epsilon/4,-yoffset,-epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    head_length = activation_head_length if indirect2 == 'increasing' else inhibition_head_length
    arrow_new(ax2,i-epsilon,-yoffset-0.52*total_length_y,epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
ax2.set_yticks([])
plt.gcf().subplots_adjust(bottom=0.3,left=0.2)
plt.savefig('total_count_of_specific_ffls_nice_N%i.pdf' % N)


#proportion, specific FFL type stratified by model - sth odd is going on, not always adding up to 1 (unknown types?)
DONT_SHOW_ZERO_FFLS_NETWORKS = True
LOG=False
cmap = matplotlib.cm.tab20c
cmap = matplotlib.cm.Paired
#cmap = matplotlib.cm.tab10
sorted_sizes = np.array(sorted(zip(nr_real_ffls,nr_specific_ffl[0]/nr_real_ffls,nr_specific_ffl[1]/nr_real_ffls,nr_specific_ffl[2]/nr_real_ffls,nr_specific_ffl[3]/nr_real_ffls,nr_specific_ffl[4]/nr_real_ffls,nr_specific_ffl[5]/nr_real_ffls,nr_specific_ffl[6]/nr_real_ffls,nr_specific_ffl[7]/nr_real_ffls,nr_unkown_ffls/nr_real_ffls),key=lambda x: (x[0],x[1]))).T
if DONT_SHOW_ZERO_FFLS_NETWORKS:
    index = list(sorted_sizes[0]>0).index(True)
    sorted_sizes = sorted_sizes[:,index:]
f,ax=plt.subplots(figsize=(10,4))
ind = np.arange(sorted_sizes.shape[1])
count_coherent=0
count_incoherent=0
for kk in range(1,8+1):
    COHERENT = int(LEGEND_COH[kk-1])
    color = cmap(count_coherent) if COHERENT else cmap(4+count_incoherent)#
    #color = [1,0+0.3*count_coherent,0+0.3*count_coherent] if COHERENT else [0+0.3*count_incoherent,0+0.3*count_incoherent,1]
    ax.bar(ind,sorted_sizes[kk],width=0.7,bottom=sum(sorted_sizes[1:kk]),color=color,label=str(kk)+', '+np.array(['INCOH','COH'])[COHERENT],log=LOG)
    if COHERENT:
        count_coherent+=1
    else:
        count_incoherent+=1
kk=9
ax.bar(ind,sorted_sizes[kk],width=0.7,bottom=sum(sorted_sizes[1:kk]),color=color_unknown)
ax.set_xlim([-1,sorted_sizes.shape[1]])
ax.set_ylabel('Proportion of specific FFLs')
ax.set_xticks([])
ax.set_xlabel('Gene regulatory networks')
ax2 = ax.twinx()
ax2.semilogy(ind,sorted_sizes[0],'x-',lw=2,color=[0.2,0.9,0.2])
ax2.set_ylabel('Number of FFLs',color=[0.2,0.9,0.2])
ax2.tick_params(axis='y', labelcolor=[0.2,0.9,0.2])
ax.legend(loc=8)
plt.savefig('proportion_of_specific_type_ffls_N%i.pdf' % N)


#proportion, specific FFL type stratified by model - nice - sth odd is going on, not always adding up to 1 (unknown types?)
DONT_SHOW_ZERO_FFLS_NETWORKS = True
SHOW_LEGEND=False
LOG=False
cmap = matplotlib.cm.Paired
color_ax2 = 'k'
sorted_sizes = np.array(sorted(zip(nr_real_ffls,nr_specific_ffl[0]/nr_real_ffls,nr_specific_ffl[1]/nr_real_ffls,nr_specific_ffl[2]/nr_real_ffls,nr_specific_ffl[3]/nr_real_ffls,nr_specific_ffl[4]/nr_real_ffls,nr_specific_ffl[5]/nr_real_ffls,nr_specific_ffl[6]/nr_real_ffls,nr_specific_ffl[7]/nr_real_ffls,nr_unkown_ffls/nr_real_ffls),key=lambda x: (x[0],x[1]))).T
if DONT_SHOW_ZERO_FFLS_NETWORKS:
    index = list(sorted_sizes[0]>0).index(True)
    sorted_sizes = sorted_sizes[:,index:]
f,ax=plt.subplots(figsize=(10,4))
ind = np.arange(sorted_sizes.shape[1])
bottom = np.zeros(sorted_sizes[order[0]].shape)
for i in range(8):
    ax.bar(ind,sorted_sizes[1+order[i]],bottom=bottom,color=cmap(i),log=LOG,width=0.7)
    bottom +=sorted_sizes[1+order[i]]     
kk=9
ax.bar(ind,sorted_sizes[kk],width=0.7,bottom=sum(sorted_sizes[1:kk]),color=color_unknown)
ax.set_xlim([-1,sorted_sizes.shape[1]])
ax.set_ylabel('Proportion of specific FFLs')
ax.set_xticks([])
ax.set_xlabel('Gene regulatory networks')
ax2 = ax.twinx()
ax2.semilogy(ind,sorted_sizes[0],'-',lw=3,color=color_ax2)
ax2.set_ylabel('Number of FFLs',color=color_ax2)
ax2.tick_params(axis='y', labelcolor=color_ax2)
if SHOW_LEGEND:
    yoffset = -1.26
    epsilon = 0.25*((sorted_sizes.shape[1]+1)/9)
    total_length_y = 0.2
    activation_head_length = 7
    inhibition_head_length = 0.1
    ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
    ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
    ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
    ax.set_ylim([0,1])
    order = np.append(np.arange(8)[LEGEND_COH],np.arange(8)[~LEGEND_COH])
    for i,ii in enumerate(np.linspace(-1,sorted_sizes.shape[1],10)[1:-1]):
        direct,indirect1,indirect2 = LEGEND[order[i]]
        head_width = 7
        color=cmap(i)
        head_length = activation_head_length if direct == 'increasing' else inhibition_head_length
        arrow_new(ax,ii+epsilon,-yoffset,0,-total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        head_length = activation_head_length if indirect1 == 'increasing' else inhibition_head_length
        arrow_new(ax,ii+epsilon/4,-yoffset,-epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
        head_length = activation_head_length if indirect2 == 'increasing' else inhibition_head_length
        arrow_new(ax,ii-epsilon,-yoffset-0.52*total_length_y,epsilon*1.25,-0.48*total_length_y,fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    plt.gcf().subplots_adjust(top=0.75)
plt.savefig('proportion_of_specific_type_ffls_nice_N%i.pdf' % N,bbox_inches = "tight")

pmids_models_loaded = [int(''.join([entry for entry in el.split('.')[0].split('_') if entry[0] in '123456789'])) for el in models_loaded]

pd.DataFrame()




## find all FFLs and determine if the direct or the indirect regulator of the common target is more important in terms of canalization
def build_dummy_vector(kis):
    res = []
    for i,el in enumerate(kis):
        res.extend([i+1]*el)
    return res

all_ffls = []
type_ffl = []
res = []
res_detailed = []
nr_layers = 5 #set this or run once and then see what the max is
res_mat = np.zeros((nr_layers+1,nr_layers+1),dtype=int)
for i in range(len(Fs)):
    F = Fs[i]
    I = Is[i]
    A = can.adjacency_matrix(I,constantss[i])
    (ffls,types) = can.get_ffls(A,F,I)
    #res.append([])
    for [regulator,intermediate,target] in ffls:
        f = F[target]
        index_regulator = list(I[target]).index(regulator)
        index_intermediate = list(I[target]).index(intermediate)
        (n,k,can_inputs,can_outputs,corefunction,can_order) = can.get_canalizing_depth_inputs_outputs_corefunction_order(f)
        kis = can.get_layer_structure_given_outputs_corefunction(can_outputs,corefunction,n)
        dummy_kis = build_dummy_vector(kis)

        try:
            index_in_canalizing_order_regulator = list(can_order).index(index_regulator)
            layer_regulator = dummy_kis[index_in_canalizing_order_regulator]
        except ValueError:
            layer_regulator = 1000
        try:
            index_in_canalizing_order_intermediate = list(can_order).index(index_intermediate)
            layer_intermediate = dummy_kis[index_in_canalizing_order_intermediate]
        except ValueError:
            layer_intermediate = 2000
        layers = [layer_regulator-1 if layer_regulator<1000 else nr_layers,layer_intermediate-1 if layer_intermediate<2000 else nr_layers]
        res_mat[layers[0],layers[1]] += 1
        res_detailed.append(','.join(list(map(str,layers))))
        res.append(layer_regulator-layer_intermediate)
res = np.array(res)
res_detailed = np.array(res_detailed)

pd.DataFrame(res_mat,columns=['Layer %i' % (i+1) for i in range(nr_layers)]+['not canalizing'],index=['Layer %i' % (i+1) for i in range(nr_layers)]+['not canalizing']).to_excel('layer_regulator_row_vs_layer_intermediate_col_in_ffls_N%i.xlsx' % len(Fs))
 
    

























## Feedbackloops and triads
max_loop=6
all_loops = []
all_types = []
triads = []
n_variables = [len(variables) for variables in variabless]
for i in range(len(Fs)):
    if all_max_degrees[i]>max_degree:
        all_types.append([])
        all_loops.append([])
        continue
    F = Fs[i]
    I = Is[i]
    degree = degrees[i]
    edges = []
    for j,regulators in enumerate(I):
        if j>=n_variables[i]: #exclude constant self-loops
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

#NR FFL vs Network Size
f,ax = plt.subplots()
for k in range(max_loop):
    ax.scatter(n_variables,nr_real_loops[k],alpha=0.5,label='%i-loops' % (k+1))
ax.legend(loc='best')

#NR FFL vs Average degree
f,ax = plt.subplots()
for k in range(max_loop):
    ax.scatter(all_avg_degrees,nr_real_loops[k],alpha=0.5,label='%i-loops' % (k+1))
ax.legend(loc='best')

#count positive vs negative vs unknown loops, stratified by model
color_unknown = [0.5,0.5,0.5]
color_neg = [0.7,0.7,1]
color_pos = [1,0.7,0.7]

for k in range(max_loop):
    DONT_SHOW_ZERO_LOOP_NETWORKS = True
    LOG=False
    sorted_sizes = np.array(sorted(zip(nr_real_loops[k],nr_pos_loops[k],nr_neg_loops[k],nr_unknown_loops[k]),key=lambda x: (x[0],x[1]))).T
    if DONT_SHOW_ZERO_LOOP_NETWORKS:
        index = list(sorted_sizes[0]>0).index(True)
        sorted_sizes = sorted_sizes[:,index:]
    f,ax=plt.subplots(figsize=(10,4))
    ind = np.arange(sorted_sizes.shape[1])
    ax.bar(ind,sorted_sizes[1],width=0.7,color=color_pos,label='positive %i-loops' % (k+1),log=LOG)
    ax.bar(ind,sorted_sizes[2],width=0.7,bottom=sorted_sizes[1],color=color_neg,label='negative %i-loops' % (k+1),log=LOG)
    ax.bar(ind,sorted_sizes[3],width=0.7,bottom=sorted_sizes[1]+sorted_sizes[2],color=color_unknown,label='unknown type %i-loops' % (k+1),log=LOG)
    ax.set_xlim([-1,sorted_sizes.shape[1]])
    ax.set_ylabel('Number of %i-loops' % (k+1))
    ax.set_xticks([])
    ax.set_xlabel('Gene regulatory networks')
    ax.legend(loc=2)
    plt.savefig('number_of_pos_neg_unknown_%iloops_N%i.pdf' % (k+1,N))

#total count positive vs negative vs unknown loops
for k in range(max_loop):
    f,ax = plt.subplots()
    ax.bar([0],[sum(nr_pos_loops[k])],color=color_pos)
    ax.bar([1],[sum(nr_neg_loops[k])],color=color_neg)
    ax.bar([2],[sum(nr_unknown_loops[k])],color=color_unknown)
    ax.set_xticks(range(3))
    ax.set_xticklabels(['positive','negative','unknown'],rotation=90)
    ax.set_ylabel('total number of %i-loops' % (k+1))
    plt.tight_layout()
    plt.savefig('total_count_of_%iloops_N%i.pdf' % (k+1,N))

#total count, specific type of loops
cmap = matplotlib.cm.Set1
for k in range(max_loop): 
    f,ax = plt.subplots()
    for i in range(k+2):
        ax.bar([i],[sum(nr_specific_k_loops[k,i,:])],color=cmap(i))
    ax.bar([k+2],[sum(nr_unknown_loops[k])],color=color_unknown)
    ax.set_xticks(range(k+3))
    ax.set_xticklabels([str(k+1-i)+r'$+$'+'\n'+str(i)+r'$-$' for i in range(k+2)]+['N.A.'])
    #ax.set_xticklabels([(r'$+$'*(k+1-i))+(r'$-$'*(i)) for i in range(k+2)]+['unknown'],rotation=90)
    ax.set_ylabel('total number of %i-loops' % (k+1))
    plt.tight_layout()
    plt.savefig('total_count_of_specific_%iloops_N%i.pdf' % (k+1,N))

#expected proportion based on proportion pos vs neg regulations
def associate_number_monotonic(*args):
    res = []
    for el in args:
        if el == 'decreasing':
            res.append(1)
        elif el== 'increasing':
            res.append(0)
        elif el== 'not essential':
            res.append(-2)        
        else:
            res.append(-1)        
    return res

res = []
res_separate = []
for i,F in enumerate(Fs):
    dummy = []
    for f in F:
        if len(f)==0: #happens if actual degree_f > max_degree
            continue
        (NONDEGENERATED,monotonic) = can.is_monotonic(f,True)
        nrs = associate_number_monotonic(*monotonic)
        n_essential=sum([el!=-2 for el in nrs])
        n = len(nrs)
        for el in nrs:
            res.append([i,n_essential,n,el])
            if el in [0,1]:
                dummy.append(el)
    res_separate.append((len(dummy)-sum(dummy))/len(dummy))

res = np.array(res)
#np.c_[np.array(np.round(res_separate,2)*100,dtype=int),nr_specific_k_loops[5,:-1,:].T]

prop_pos = sum(res[:,3]==0)/(sum(res[:,3]==0)+sum(res[:,3]==1))
neg_of_ns = [[4,4],[3,4],[6,6],[5,6],[4,6]]

f,ax = plt.subplots(figsize=(9,4))
width=0.8
height_rectangle = 0.08
epsilon = 0.05
colors = ['blue','orange']
for i,neg_of_n in enumerate(neg_of_ns):
    nr_neg,k = neg_of_n
    nr_pos = k-nr_neg
    more_pos = sum(nr_specific_k_loops[k-1,nr_pos,:])
    more_neg = sum(nr_specific_k_loops[k-1,nr_neg,:])
    total_nr_loops = more_pos+more_neg
    ax.bar([2*i-width/2-epsilon],[more_pos/total_nr_loops],width=width,color=colors[0],alpha=0.3)
    ax.bar([2*i-width/2-epsilon],[more_neg/total_nr_loops],bottom=[more_pos/total_nr_loops],width=width,color=colors[1],alpha=0.3)
    ax.text(2*i,1.05,'n=%i' % total_nr_loops,va='center',ha='center',clip_on=False )
    ax.text(2*i-width/2-epsilon,more_pos/total_nr_loops/2,str(nr_neg)+r'$+$'+'\n'+str(nr_pos)+r'$-$',va='center',ha='center')
    ax.text(2*i-width/2-epsilon,more_pos/total_nr_loops+more_neg/total_nr_loops/2,str(nr_pos)+r'$+$'+'\n'+str(nr_neg)+r'$-$',va='center',ha='center')
    expected_more_pos = prop_pos**nr_neg*(1-prop_pos)**nr_pos
    expected_more_neg = prop_pos**nr_pos*(1-prop_pos)**nr_neg
    total_expected = expected_more_pos+expected_more_neg
    ax.bar([2*i+width/2+epsilon],[expected_more_pos/total_expected],width=width,color=colors[0])
    ax.bar([2*i+width/2+epsilon],[expected_more_neg/total_expected],bottom=[expected_more_pos/total_expected],width=width,color=colors[1])

ax.add_patch(matplotlib.patches.Rectangle([2*len(neg_of_ns),0.75],1,height_rectangle,color=colors[0],clip_on=False))
ax.add_patch(matplotlib.patches.Rectangle([2*len(neg_of_ns),0.65],1,height_rectangle,color=colors[1],clip_on=False))
ax.text(2*len(neg_of_ns)+1.2,0.75+height_rectangle/2,'more positive regulations',ha='left',va='center')
ax.text(2*len(neg_of_ns)+1.2,0.65+height_rectangle/2,'more negative regulations',ha='left',va='center')

ax.add_patch(matplotlib.patches.Rectangle([2*len(neg_of_ns),0.35],0.5,height_rectangle,color=colors[0],alpha=0.3,clip_on=False))
ax.add_patch(matplotlib.patches.Rectangle([2*len(neg_of_ns)+0.5,0.35],0.5,height_rectangle,color=colors[1],alpha=0.3,clip_on=False))
ax.add_patch(matplotlib.patches.Rectangle([2*len(neg_of_ns),0.25],1,height_rectangle,color=colors[0],alpha=1,clip_on=False))
ax.add_patch(matplotlib.patches.Rectangle([2*len(neg_of_ns)+0.5,0.25],0.5,height_rectangle,color=colors[1],alpha=1,clip_on=False))
ax.text(2*len(neg_of_ns)+1.2,0.35+height_rectangle/2,'observed',ha='left',va='center')
ax.text(2*len(neg_of_ns)+1.2,0.25+height_rectangle/2,'expected',ha='left',va='center')

#ax.set_ylim([0,1])
## Hide the right and top spines
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
## Only show ticks on the left and bottom spines
#ax.yaxis.set_ticks_position('left')
#ax.xaxis.set_ticks_position('none')
#ax.xaxis.set_ticks([])
#ax.set_ylabel('Proportion')
#plt.savefig('pos_vs_neg_regulations_in_loop_N%i_nice.pdf' % len(Fs),bbox_inches = "tight")

arrow_new(ax,2*len(neg_of_ns)+1,1.05,-1,0,clip_on=False)
ax.text(2*len(neg_of_ns)+1.2,1.05,'total number observed',ha='left',va='center')

arrow_new(ax,2*len(neg_of_ns)+1,-0.07,-1,0,clip_on=False)
ax.text(2*len(neg_of_ns)+1.2,-0.07,'type of feedback loop',ha='left',va='center')

ax.plot([-width-epsilon,2*1+width+epsilon],[-0.13,-0.13],'k-',clip_on=False,lw=0.5)
ax.text(1,-0.2,'4-loops',ha='center',va='center')
ax.plot([2*2-width-epsilon,2*4+width+epsilon],[-0.13,-0.13],'k-',clip_on=False,lw=0.5)
ax.text(6,-0.2,'6-loops',ha='center',va='center')

ax.set_ylim([0,1])
ax.set_xlim([-1,2*len(neg_of_ns)-1])
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('none')
ax.xaxis.set_ticks(2*np.arange(len(neg_of_ns)))
ax.xaxis.set_ticklabels(['positive' if el[0]%2==0 else 'negative' for el in neg_of_ns])
ax.set_ylabel('Proportion')
plt.savefig('pos_vs_neg_regulations_in_loop_N%i_nice.pdf' % len(Fs),bbox_inches = "tight")














## Clustering of two FFLs or a FFL and a 2-loop or a FFL and a 3-loop
all_ffls = []
type_ffl_ffls = []
type_ffl_2l = []
type_ffl_3l = []

DEBUG = True

for ii in range(N):
    if all_max_degrees[ii]>max_degree:
        type_ffl_ffls.append([])
        continue
    F = Fs[ii]
    if len(F)>400:
        continue
    I = Is[ii]
    degree = degrees[ii]
    n_variables = len(variabless[ii])
    constants = constantss[ii]
    A = can.adjacency_matrix(Is[ii],constantss[ii])
    (ffls,types) = can.get_ffls(A,F,I)
    
    type_ffl_ffl = []
    all_ffls.append(list(map(can.get_ffl_type_number,types)))
    #Determines types of motifs from 2018 Science Advances paper and new motifs (FFL and 2-loop)
    n_ffls = len(ffls)
    for i in range(n_ffls):
        for j in range(i+1,n_ffls):
            intersect = set.intersection(set(ffls[i]),set(ffls[j]))
            if len( intersect )>0:
                if len(intersect)==3: #new types not included in 2018 paper
                    types = dict(zip(ffls[i],['r','i','t']))
                    for new_type,el in zip(['r','i','t'],ffls[j]):
                        try:
                            type_el = types[el]
                            if new_type != type_el:
                                types.update({el:'m'})
                        except:
                            types.update({el:new_type})
                    sorted_types = list(types.values())
                    sorted_types.sort()
                    string_type = ''.join(sorted_types)
                    if DEBUG:
                        print(string_type)
                    if string_type == 'mmr':
                        type_ffl_ffl.append(13)
                    elif string_type == 'mmt':
                        type_ffl_ffl.append(14)
                    elif string_type == 'mmm': #When all the nodes have switched roles
                        if ffls[i][1] == ffls[j][0] or ffls[i][1] == ffls[j][2]:
                            type_ffl_ffl.append(15)
                        else:
                            type_ffl_ffl.append(17)
                elif len(intersect)==2: #types 7-12 from 2018 paper
                    types = dict(zip(ffls[i],['r','i','t']))
                    for new_type,el in zip(['r','i','t'],ffls[j]):
                        try:
                            type_el = types[el]
                            if new_type != type_el:
                                types.update({el:'m'})
                        except:
                            types.update({el:new_type})
                    sorted_types = list(types.values())
                    sorted_types.sort()
                    string_type = ''.join(sorted_types)
                    if DEBUG:
                        print(string_type)
                    if string_type == 'irrt':
                        type_ffl_ffl.append(7)
                    elif string_type == 'imrt':
                        if ffls[i][2] == ffls[j][2]:#if the target in both FFLs is the same
                            type_ffl_ffl.append(8)
                        else:
                            type_ffl_ffl.append(11)
                    elif string_type == 'mmrt':
                        type_ffl_ffl.append(9)
                    elif string_type == 'iirt':
                        type_ffl_ffl.append(10)                
                    elif string_type == 'irtt':
                        type_ffl_ffl.append(12)
                elif len(intersect)==1: #types 1-6 from 2018 paper
                    types = dict(zip(ffls[i],['r','i','t']))
                    for new_type,el in zip(['r','i','t'],ffls[j]):
                        try:
                            type_el = types[el]
                            if new_type != type_el:
                                types.update({el:'m'})
                        except:
                            types.update({el:new_type})
                    sorted_types = list(types.values())
                    sorted_types.sort()
                    string_type = ''.join(sorted_types)
                    if DEBUG:
                        print(string_type)
                    if string_type == 'iirrt':
                        type_ffl_ffl.append(1)
                    elif string_type == 'imrrt':
                        type_ffl_ffl.append(2)
                    elif string_type == 'iimrt':
                        type_ffl_ffl.append(3)
                    elif string_type == 'irrtt':
                        type_ffl_ffl.append(4)                
                    elif string_type == 'imrtt':
                        type_ffl_ffl.append(5)
                    elif string_type == 'iirtt':
                        type_ffl_ffl.append(6)                        
    type_ffl_ffls.append(type_ffl_ffl)
    
    #FFLs clustered with 3 loops
    max_loop_length = 3
    G=can.generate_networkx_graph_from_edges(I,n_variables)
    loops = list(can.simple_cycles(G,max_len = max_loop_length))
    two_loops = [el for el in loops if len(el)==2]
    three_loops = [el for el in loops if len(el)==3]
    
    for i in range(n_ffls):
        for j in range(len(two_loops)):            
            intersect = set.intersection(set(ffls[i]),set(two_loops[j]))
            if len( intersect )>0:
                types = dict(zip(ffls[i],['r','i','t']))
                for el in list(intersect):
                    types.update({el:'m'})
                sorted_types = list(types.values())
                sorted_types.sort()
                string_type = ''.join(sorted_types) 
                print(string_type)
                if string_type == 'mmr':
                    type_ffl_2l.append(13)
                elif string_type == 'mmt':
                    type_ffl_2l.append(14)
                elif string_type == 'imm':
                    type_ffl_2l.append(18)
                elif string_type == 'imt':
                    type_ffl_2l.append(19)
                elif string_type == 'mrt':
                    type_ffl_2l.append(20)
                elif string_type == 'imr':
                    type_ffl_2l.append(21)

all_type_ffl_ffls = []
for i,el in enumerate(type_ffl_ffls): 
    all_type_ffl_ffls.extend(el)
    
all_types = list(set(all_type_ffl_ffls))
n_types = len(all_types)
dict_types = dict(zip(all_types,range(n_types)))

nr_ffls = np.array([len(el) for el in all_ffls])
nr_specific_ffl_ffl = np.zeros((n_types,N),dtype=int)
for i,el in enumerate(type_ffl_ffls): 
    for typ in el:
        nr_specific_ffl_ffl[dict_types[typ],i]+=1
nr_total_ffl_ffl = np.sum(nr_specific_ffl_ffl,0)
     
sum(np.sum(nr_specific_ffl_ffl,1)[:6])
sum(np.sum(nr_specific_ffl_ffl,1)[6:12])
sum(np.sum(nr_specific_ffl_ffl,1)[12:])


print(pd.value_counts(total_type_ffl_ffl))
print(pd.value_counts(total_type_ffl_2l))


#total count, specific FFL type - nice
import matplotlib
def arrow_new(self, x, y, dx, dy, **kwargs):
    kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
    kwargs.setdefault('fc', 'black')
    x = self.convert_xunits(x)
    y = self.convert_yunits(y)
    dx = self.convert_xunits(dx)
    dy = self.convert_yunits(dy)
    posA = x, y
    posB = x+dx, y+dy
    a = matplotlib.patches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
    self.add_artist(a)
    return a

types_graphic = [list([[0,1],[0,2],[1,2],[3,2],[3,4],[4,2]]),
                 list([[0,1],[0,2],[2,1],[3,2],[3,4],[4,2]]),
                 list([[0,1],[0,2],[1,2],[2,3],[3,4],[2,4]]),
                 list([[0,1],[0,2],[2,1],[3,2],[3,4],[2,4]]),
                 list([[0,1],[2,0],[2,1],[3,2],[3,4],[2,4]]),
                 list([[0,1],[2,0],[2,1],[2,3],[3,4],[2,4]]),
                 list([[0,1],[0,2],[1,2],[3,1],[3,2]]),
                 list([[0,1],[0,2],[1,2],[1,3],[3,2]]),
                 list([[0,1],[0,2],[1,2],[1,3],[2,3]]),
                 list([[1,0],[0,2],[1,2],[1,3],[3,2]]),
                 list([[1,0],[2,0],[1,2],[1,3],[3,2]]),
                 list([[1,0],[2,0],[1,2],[1,3],[2,3]]),
                 list([[0,1],[0,2],[1,2],[2,1]]),
                 list([[0,1],[0,2],[1,2],[1,0]]),
                 list([[0,1],[0,2],[1,2],[1,0],[2,0]])]
type_graphic = types_graphic[0]

n_types=15
order=list(range(n_types))

cmap = matplotlib.cm.tab20

COLOR_BY_TYPE = False
colors_r_i_t_m = [[0.3,0.3,1],[1,0.9,0.2],[0.8,0,0],[0.5,0.5,0.5]]
f,ax = plt.subplots(figsize=(11,4))
width=0.8
for i in range(n_types):
    ax.bar([i],[sum(nr_specific_ffl_ffl[order[i],:])],color=cmap(i),width=width)
ax.set_xticks(range(n_types))
ax.set_xticklabels(['' for i in range(n_types)],rotation=90)
ax.xaxis.set_ticks_position('none') 
ax.set_ylabel('total count')
ax.set_xlim([-.5,n_types-.5])

ax2 = ax.twinx()

ax2.set_ylim([0,1])    
delta_y = 0.13
delta_x = 0.4*width
total_length_y = 0.2
activation_head_length = 7
inhibition_head_length = 0.1
#ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
#ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
#ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
#ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
#ax.set_ylim([0,1])

head_width = 4
head_length = 4

ycenter = -0.17
points = [(-1,1),(1,1),(0,0),(-1,-1),(1,-1)] #from top to bottom, left to right
for j in range(6):
    xcenter=j
    color=cmap(j) if COLOR_BY_TYPE else 'k'
    type_graphic = types_graphic[j]
    for i in range(len(type_graphic)):
        #ax.text(xcenter-0.8*delta_x,ycenter,str(all_types[j]),color=color,va='center',ha='center')
        arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    for ii,point in enumerate(points):
        if not COLOR_BY_TYPE:
            ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
            if ii not in ffls_in_cluster[:,1:]: #regulator
                color=colors_r_i_t_m[0]
            elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                color=colors_r_i_t_m[1]
            elif ii not in ffls_in_cluster[:,:2]: #target
                color=colors_r_i_t_m[2]
            else:
                color=colors_r_i_t_m[3] 
        ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=6,clip_on=False)
points = [(0,1),(-1,0),(1,0),(0,-1)] #from top to bottom, left to right
for j in range(6,12):
    xcenter=j
    color=cmap(j) if COLOR_BY_TYPE else 'k'
    type_graphic = types_graphic[j]
    for i in range(len(type_graphic)):
        #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
        arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    for ii,point in enumerate(points):
        if not COLOR_BY_TYPE:
            ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
            if ii not in ffls_in_cluster[:,1:]: #regulator
                color=colors_r_i_t_m[0]
            elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                color=colors_r_i_t_m[1]
            elif ii not in ffls_in_cluster[:,:2]: #target
                color=colors_r_i_t_m[2]
            else:
                color=colors_r_i_t_m[3]                
        ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=6,clip_on=False)
points = [(0,.5),(-1,-.5),(1,-.5)] #from top to bottom, left to right
for j in range(12,15):
    xcenter=j
    color=cmap(j) if COLOR_BY_TYPE else 'k'
    type_graphic = types_graphic[j]
    for i in range(len(type_graphic)):
        #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
        arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    for ii,point in enumerate(points):
        if not COLOR_BY_TYPE:
            ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
            if ii not in ffls_in_cluster[:,1:]: #regulator
                color=colors_r_i_t_m[0]
            elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                color=colors_r_i_t_m[1]
            elif ii not in ffls_in_cluster[:,:2]: #target
                color=colors_r_i_t_m[2]
            else:
                color=colors_r_i_t_m[3] 
        ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=6,clip_on=False)
ax2.set_yticks([])
plt.gcf().subplots_adjust(bottom=0.3,left=0.12,right=0.99)
plt.savefig('total_count_of_specific_ffls_ffls_nice_N%i.pdf' % (N))


#nice for publication
f,ax = plt.subplots(figsize=(11,4.5))
width=0.8
for i in range(n_types):
    ax.bar([i],[sum(nr_specific_ffl_ffl[order[i],:])],color=cmap(i),width=width)
ax.set_xticks(range(n_types))
ax.set_xticklabels(['' for i in range(n_types)],rotation=90)
ax.xaxis.set_ticks_position('none') 
ax.set_ylabel('total count')
ax.set_xlim([-.5,n_types-.5])

ax2 = ax.twinx()

ax2.set_ylim([0,1])    
delta_y = 0.13
delta_x = 0.4*width
total_length_y = 0.2
activation_head_length = 7
inhibition_head_length = 0.1
#ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
#ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
#ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
#ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
#ax.set_ylim([0,1])

head_width = 4
head_length = 4

ycenter = -0.31
points = [(-1,1),(1,1),(0,0),(-1,-1),(1,-1)] #from top to bottom, left to right
for j in range(15):
    ax2.text(j,-0.1,str(j+1),va='center',ha='center',color=cmap(j),clip_on=False)
for j in range(6):
    xcenter=j
    color=cmap(j) if COLOR_BY_TYPE else 'k'
    type_graphic = types_graphic[j]
    for i in range(len(type_graphic)):
        #ax.text(xcenter-0.8*delta_x,ycenter,str(all_types[j]),color=color,va='center',ha='center')
        arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    for ii,point in enumerate(points):
        if not COLOR_BY_TYPE:
            ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
            if ii not in ffls_in_cluster[:,1:]: #regulator
                color=colors_r_i_t_m[0]
            elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                color=colors_r_i_t_m[1]
            elif ii not in ffls_in_cluster[:,:2]: #target
                color=colors_r_i_t_m[2]
            else:
                color=colors_r_i_t_m[3] 
        ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=6,clip_on=False)
points = [(0,1),(-1,0),(1,0),(0,-1)] #from top to bottom, left to right
for j in range(6,12):
    xcenter=j
    color=cmap(j) if COLOR_BY_TYPE else 'k'
    type_graphic = types_graphic[j]
    for i in range(len(type_graphic)):
        #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
        arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    for ii,point in enumerate(points):
        if not COLOR_BY_TYPE:
            ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
            if ii not in ffls_in_cluster[:,1:]: #regulator
                color=colors_r_i_t_m[0]
            elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                color=colors_r_i_t_m[1]
            elif ii not in ffls_in_cluster[:,:2]: #target
                color=colors_r_i_t_m[2]
            else:
                color=colors_r_i_t_m[3]                
        ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=6,clip_on=False)
points = [(0,.5),(-1,-.5),(1,-.5)] #from top to bottom, left to right
for j in range(12,15):
    xcenter=j
    color=cmap(j) if COLOR_BY_TYPE else 'k'
    type_graphic = types_graphic[j]
    for i in range(len(type_graphic)):
        #ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
        arrow_new(ax2,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    for ii,point in enumerate(points):
        if not COLOR_BY_TYPE:
            ffls_in_cluster = np.array(can.get_ffls_from_I(can.edgelist_to_I(type_graphic)[0]))
            if ii not in ffls_in_cluster[:,1:]: #regulator
                color=colors_r_i_t_m[0]
            elif ii not in ffls_in_cluster[:,0] and ii not in ffls_in_cluster[:,2]: #intermediate
                color=colors_r_i_t_m[1]
            elif ii not in ffls_in_cluster[:,:2]: #target
                color=colors_r_i_t_m[2]
            else:
                color=colors_r_i_t_m[3] 
        ax2.plot([xcenter+delta_x*point[0]],[ycenter+delta_y*point[1]],'o',color=color,markersize=6,clip_on=False)
ax2.set_yticks([])
plt.gcf().subplots_adjust(bottom=0.3,left=0.12,right=0.99)
plt.savefig('total_count_of_specific_ffls_ffls_nice_N%i.pdf' % (N),bbox_inches = "tight")












SHOW_LEGEND = False
DONT_SHOW_ZERO_FFLS_NETWORKS = True
LOG=False
cmap = matplotlib.cm.tab20
color_ax2 = 'k'
sorted_sizes = np.array(sorted(zip(nr_total_ffl_ffl,*nr_specific_ffl_ffl/nr_total_ffl_ffl),key=lambda x: (x[0],x[1]))).T
if DONT_SHOW_ZERO_FFLS_NETWORKS:
    index = list(sorted_sizes[0]>0).index(True)
    sorted_sizes = sorted_sizes[:,index:]
f,ax=plt.subplots(figsize=(10,6) if SHOW_LEGEND else (11,4))
ind = np.arange(sorted_sizes.shape[1])
bottom = np.zeros(sorted_sizes[order[0]].shape)
for i in range(n_types):
    ax.bar(ind,sorted_sizes[1+order[i]],bottom=bottom,color=cmap(i),log=LOG,width=0.7)
    bottom +=sorted_sizes[1+order[i]]        
ax.set_xlim([-1,sorted_sizes.shape[1]])
ax.set_ylabel('Proportion of specific FFL cluster')
ax.set_xticks([])
ax.set_xlabel('Gene regulatory networks')
ax2 = ax.twinx()
ax2.semilogy(ind,sorted_sizes[0],'-',lw=3,color=color_ax2)
ax2.set_ylabel('Number of FFL clusters',color=color_ax2)
ax2.tick_params(axis='y', labelcolor=color_ax2)
delta_y = 0.13
delta_x = 0.4*((sorted_sizes.shape[1]+1)/7)
total_length_y = 0.2
activation_head_length = 7
inhibition_head_length = 0.1
#ax.plot([-1+0.75*(sorted_sizes.shape[1]+1)/9,-1+4.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
#ax.plot([-1+4.75*(sorted_sizes.shape[1]+1)/9,-1+8.25*(sorted_sizes.shape[1]+1)/9],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
#ax.text(-1+2.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
#ax.text(-1+6.5*(sorted_sizes.shape[1]+1)/9,1.25*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
ax.set_ylim([0,1])

if SHOW_LEGEND:
    
    ycenter = 1 + 0.05 + delta_y
    points = [(-1,1),(1,1),(0,0),(-1,-1),(1,-1)] #from top to bottom, left to right
    for j,xcenter in enumerate(np.linspace(-1,sorted_sizes.shape[1],8)[1:-1]):
        head_width = 7
        head_length = 7
        color=cmap(j)
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            ax.text(xcenter-0.8*delta_x,ycenter,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    
    ycenter = 1 + 0.05 + 3.4*delta_y
    points = [(0,1),(-1,0),(1,0),(0,-1)] #from top to bottom, left to right
    for j,xcenter in enumerate(np.linspace(-1,sorted_sizes.shape[1],8)[1:-1]):
        j=j+6
        head_width = 7
        head_length = 7
        color=cmap(j)
        print(j,color)
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    
    ycenter = 1 + 0.05 + 4.8*delta_y
    points = [(0,1),(-1,0),(1,0)] #from top to bottom, left to right
    xcenters = np.linspace(-1,sorted_sizes.shape[1],8)[1:-1]
    xcenters = xcenters[0::2]+(xcenters[1]-xcenters[0])/2
    for j,xcenter in enumerate(xcenters):
        j=j+12
        head_width = 7
        head_length = 7
        color=cmap(j)
        print(j,color)
        type_graphic = types_graphic[j]
        for i in range(len(type_graphic)):
            ax.text(xcenter-0.8*delta_x,ycenter+0.8*delta_y,str(all_types[j]),color=color,va='center',ha='center')
            arrow_new(ax,xcenter+delta_x*points[type_graphic[i][0]][0]+0.1*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),ycenter+delta_y*points[type_graphic[i][0]][1]+0.1*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),0.8*delta_x*(points[type_graphic[i][1]][0]-points[type_graphic[i][0]][0]),0.8*delta_y*(points[type_graphic[i][1]][1]-points[type_graphic[i][0]][1]),fc=color,ec=color,color=color,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
                    
    plt.gcf().subplots_adjust(top=0.57,bottom=0.05)
plt.savefig('proportion_of_specific_type_ffls_clusters_nice_N%i.pdf' % (N),bbox_inches = "tight")




#Try again: 


























## FFL in complete E.Coli K-12 network
file = 'generegulation_tmp.txt'
f = open(file,'r')
text = f.read().splitlines()
f.close()
genes = []
dict_genes = dict()
count=0
I=[]
I_type = []
for el in text:
    elsplit = el.split('\t')
    if len(elsplit)==9:
        tf = elsplit[1]
        target = elsplit[7]
        try:
            tf_index = dict_genes[tf]
        except KeyError:
            I.append([])
            I_type.append([])
            tf_index = count
            dict_genes.update({tf:count})
            count+=1
        try:
            target_index = dict_genes[target]
        except KeyError:
            I.append([])
            I_type.append([])
            target_index = count
            dict_genes.update({target:count})
            count+=1
        I[target_index].append(tf_index)
        if elsplit[8] not in ['activator','repressor']:
            print(tf,target,elsplit[8])
        I_type[target_index].append('increasing' if elsplit[8]=='activator' else ('decreasing' if elsplit[8]=='repressor' else ('not monotonic' if elsplit[8]=='dual' else 'unknown')))

ffls,types = can.get_ffls_from_I(I,I_type)
types_converted = list(map(can.get_ffl_type_number,types))

LEGEND = list(itertools.product(['decreasing', 'increasing'], repeat=3))
LEGEND_NR = list(map(can.get_ffl_type_number,LEGEND))
LEGEND_COH = np.array(list(map(can.is_ffl_coherent,LEGEND)))

nr_coh_ffls = 0
nr_incoh_ffls = 0
nr_unkown_ffls = 0
nr_notreal_ffls = 0
nr_specific_ffl = np.zeros(8,dtype=int)
for el in types_converted:
    if el==-1:
        nr_unkown_ffls+=1
    elif el==-2:
        nr_notreal_ffls+=1
    elif LEGEND_COH[el] == True:
        nr_coh_ffls += 1
        nr_specific_ffl[el] += 1
    else:
        nr_incoh_ffls += 1
        nr_specific_ffl[el] += 1
        
#total count, specific FFL type - nice
order = np.append(np.arange(8)[LEGEND_COH],np.arange(8)[~LEGEND_COH])
cmap = matplotlib.cm.Paired

f,ax = plt.subplots()
width=0.8
for i in range(8):
    ax.bar([i],[(nr_specific_ffl[order[i]])],color=cmap(i),width=width)
ax.bar([8],[(nr_unkown_ffls)],color=color_unknown)
ax.set_xticks(range(9))
ax.set_xticklabels(['' for i in range(8)]+['unknown'],rotation=90)
ax.xaxis.set_ticks_position('none') 
ax.set_ylabel('total count')
ax2 = ax.twinx()
yoffset = 0.05
epsilon = 0.25
total_length_y = 0.2
activation_head_length = 7
inhibition_head_length = 0.1
ax2.plot([0-width/2,3+width/2],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
ax2.plot([4-width/2,7+width/2],[1.2*(-yoffset-total_length_y),1.2*(-yoffset-total_length_y)],'k-',clip_on=False)
ax2.text(1.5,1.45*(-yoffset-total_length_y),'coherent FFLs',va='center',ha='center',clip_on=False)
ax2.text(5.5,1.45*(-yoffset-total_length_y),'incoherent FFLs',va='center',ha='center',clip_on=False)
ax2.set_ylim([0,1])
for i in range(8):
    direct,indirect1,indirect2 = LEGEND[order[i]]
    head_width = 7
    head_length = activation_head_length if direct == 'increasing' else inhibition_head_length
    arrow_new(ax2,i+epsilon,-yoffset,0,-total_length_y,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Simple(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    head_length = activation_head_length if indirect1 == 'increasing' else inhibition_head_length
    arrow_new(ax2,i+epsilon/4,-yoffset,-epsilon*1.25,-0.48*total_length_y,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
    head_length = activation_head_length if indirect2 == 'increasing' else inhibition_head_length
    arrow_new(ax2,i-epsilon,-yoffset-0.52*total_length_y,epsilon*1.25,-0.48*total_length_y,clip_on=False,arrowstyle=matplotlib.patches.ArrowStyle.Fancy(head_length=head_length, head_width=head_width, tail_width=.4))#),head_width=0.2,head_length=head_length,length_includes_head=True,color='k')
ax2.set_yticks([])
plt.gcf().subplots_adjust(bottom=0.3,left=0.2)
plt.savefig('total_count_of_specific_ffls_nice_regulon_db_Ecoli.pdf')

## activation vs inhibition in Ecoli k-12
I_type_flat = []
for el in I_type:
    I_type_flat.extend(el)
print(pd.value_counts(I_type_flat))