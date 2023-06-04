#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:01:32 2019

@author: ckadelka
"""

##Imports

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import canalizing_function_toolbox_v1_9 as can
import load_database11 as db

plt.rcParams.update({'font.size': 16})

## load the database, choose low max_n for quick results and to only look at small models
folders = ['update_rules_cell_collective/']#, 'update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 16
max_n=1000
[Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded] = db.load_database(folders,max_degree=max_degree,max_n=max_n)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)



pmids_models_loaded = [int(''.join([entry for entry in el.split('.')[0].split('_') if entry[0] in '123456789'])) for el in models_loaded]




#interesting observations
#1) degree vs essential degree

total_genes = sum(map(len,variabless))
print( total_genes )
print( sum(map(len,constantss)) )

print( np.mean(list(map(len,variabless))) )
print( np.median(list(map(len,variabless))) )
print(sum(np.array(list(map(len,constantss)))>0))
print(max(map(len,variabless)))
print(min(map(len,variabless)))

n_variables = np.array([len(variables) for variables in variabless])
total_functions = 0
n_constants = np.array([len(constants) for constants in constantss])
vec_mean = []
vec_total = []
vec_max = []
for i in range(N):
    a=np.array(degrees[i][:n_variables[i]])
    b=np.array(degrees_essential[i][:n_variables[i]])
    dummy = np.array([np.array(el1)-np.array(el2) for el1,el2 in zip(a[a<=max_degree],b[a<=max_degree])])
    total_functions += len(dummy)
    vec_mean.append(np.mean(dummy>0))
    vec_total.append(np.sum(dummy>0))
    vec_max.append(np.max(dummy))
argmax = np.argmax(vec_max)
argmax_of_this = np.argmax(np.array(degrees[argmax])-np.array(degrees_essential[argmax]))
Fs[argmax][argmax_of_this]
#'rule 6 (i.e 7) in IGVH mutations in chronic lymphocytic leukemia_26088082.txt'
#Look at BMI1
#PTEN = ( EGR1 & FRK & INPP5D & RRM1  & ! TGFBR3 & BMI1 & ! AEBP1 & ! BMI1 & ILK)| (!ADM & !CSNK2A2  & EGR1 & FRK & INPP5D & RRM1 & TNF & ! TGFBR3 & IFNGR2 & BMI1 & !AEBP1 & ! BMI1 & ILK)

print( sum([el>0 for el in vec_total]))
print(sum(vec_total),'/',total_functions,': ',sum(vec_total)/total_functions)


f,ax=plt.subplots()
for i in range(N):
    a=np.array(degrees[i][:n_variables[i]])
    b=np.array(degrees_essential[i][:n_variables[i]])
    ax.plot(a[a<=max_degree],b[a<=max_degree],'ko',alpha=0.5)
ax.set_xlabel('Number of regulators')
ax.set_ylabel('Number of essential regulators')
ax.set_xlim([ax.get_ylim()[0],ax.get_xlim()[1]])
plt.gcf().subplots_adjust(bottom=0.2)
plt.savefig('degree_vs_actual_degree_N%i.pdf' % len(Fs))    


avg_connectivity = [np.mean(np.array(el)[np.array(el)>0]) for el in degrees_essential]

x=np.array(n_variables)
y = np.array(avg_connectivity)
from scipy import polyfit
import scipy.stats as stats
#b, m = polyfit(np.log(x), y, 1)
m,b,r,p,stderr = stats.linregress(np.log(x), y)
r,p = stats.spearmanr(x, y)
f,ax=plt.subplots()
ax.semilogx(n_variables,avg_connectivity,'ko',markersize=8,alpha=0.5)
ax.set_xlabel('Number of genes')
ax.set_ylabel('Average connectivity')
#ax.semilogx(x, b + m * np.log(x), 'x')
(x1,x2) = ax.get_xlim()
(y1,y2) = ax.get_ylim()
#ax.semilogx([x1,x2],b+m*np.log([x1,x2]),'--',color='r',lw=2)
ax.set_xlim([x1,x2])
ax.text(np.exp(np.log(x1)+0.03*(np.log(x2)-np.log(x1))),y1+0.87*(y2-y1),(r'r$=%s$' % str(np.round(r,2))) +'\n'+('p$=%s$' % str(np.round(p,2))),color='red',va='center',ha='left')
#print(stats.pearsonr(x,y))
plt.gcf().subplots_adjust(bottom=0.2,left=0.2)
plt.savefig('size_vs_avg_connectivity_N%i.pdf' % N,bbox_inches = "tight")



#models of interest for discussion with brandy et al
import itertools
def bool_to_poly_v2(f,I,variables,constants,x=[]):
    variables_and_constants  = list(variables)+list(constants)
    len_f = len(f)
    n=int(np.log2(len_f))
    if x==[]: #to reduce run time, this should be calculated once and then passed as argument
        x = list(itertools.product([0, 1], repeat = n))
    num_values = 2**n   
    text = []
    for i in range(num_values):
        if f[i]==True:
            #monomial = '*'.join([('x%i' % (j+1)) if entry==1 else ('(1-x%i)' % (j+1)) for j,entry in enumerate(x[i])])
            monomial = '*'.join([variables_and_constants[I[j]] if entry==1 else ('(1-%s)' % variables_and_constants[I[j]]) for j,entry in enumerate(x[i])])
            text.append(monomial)
    if text!=[]:
        return ' + '.join(text)
    else:
        return '0'

which = np.bitwise_and(np.array(n_variables)<25,np.array(avg_connectivity)>3)
ids = np.array(pmids_models_loaded)[which]
for i in np.arange(len(n_variables))[which]:
    F_poly = []
    F_poly_v2 = []
    for j,(f,var) in enumerate(zip(Fs[i],variabless[i])):
        F_poly.append(var + '=' + can.bool_to_poly(f))
        F_poly_v2.append(var + '=' + bool_to_poly_v2(f,Is[i][j],variabless[i],constantss[i]))
    g = open(models_loaded[i].replace('.txt','_poly.txt'),'w')
    g.write('\n'.join(F_poly))
    g.close()
    g = open(models_loaded[i].replace('.txt','_poly_withnames.txt'),'w')
    g.write('\n'.join(F_poly_v2))
    g.close()



f,ax=plt.subplots()
ax.loglog(n_variables,n_constants,'ko',markersize=8,alpha=0.5)
#m,b,r,p,stderr = stats.linregress(n_variables, n_constants)
(x1,x2) = ax.get_xlim()
(y1,y2) = ax.get_ylim()
#ax.semilogx([x1,x2],b+m*np.array([x1,x2]),'--',color='r',lw=2)
r,p = stats.spearmanr(n_variables,n_constants)
ax.set_xlim([x1,x2])
ax.text(np.exp(np.log(x1)+0.03*(np.log(x2)-np.log(x1))),y1+0.87*(y2-y1),(r'r$=%s$' % str(np.round(r,2))) +'\n'+('p$=%s$' % str(np.round(p,2) if p>0.01 else '{:.0e}'.format(p))),color='red',va='center',ha='left')
ax.set_xlabel('Number of genes')
ax.set_ylabel('Number of external parameters')
plt.gcf().subplots_adjust(bottom=0.2,left=0.2)
plt.savefig('number_of_genes_vs_constants_N%i.pdf' % N,bbox_inches = "tight")



f,ax=plt.subplots()
dummy_n_constants = np.array(n_constants,dtype=float)
dummy_n_constants[n_constants==0] = 0.5
ax.loglog(n_variables,dummy_n_constants,'ko',markersize=8,alpha=0.5)
#m,b,r,p,stderr = stats.linregress(n_variables, n_constants)
(x1,x2) = ax.get_xlim()
(y1,y2) = ax.get_ylim()
#ax.semilogx([x1,x2],b+m*np.array([x1,x2]),'--',color='r',lw=2)
r,p = stats.spearmanr(n_variables,n_constants)
ax.set_xlim([x1,x2])
ax.text(np.exp(np.log(x1)+0.03*(np.log(x2)-np.log(x1))),np.exp(np.log(y1)+0.87*(np.log(y2)-np.log(y1))),(r'r$=%s$' % str(np.round(r,2))) +'\n'+('p$=%s$' % str(np.round(p,2) if p>0.01 else '{:.0e}'.format(p))),color='red',va='center',ha='left')
ax.set_xlabel('Number of genes')
ax.set_ylabel('Number of external parameters')
ax.set_yticks([0.5,1,10,100])
ax.set_yticklabels(list(map(str,[0,1,10,100])))
ax.spines['left'].set_visible(False)
ax.plot([x1,x1],[y1,0.5],'k',lw=1)
ax.plot([x1,x1],[1,y2],'k',lw=1)

minoryticks = np.array(ax.get_yticks(minor=True))
minoryticks = minoryticks[minoryticks>1]
minoryticks = minoryticks[minoryticks<y2]
ax.set_yticks(minoryticks,minor=True)
#kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
#d = .02
#ax.plot((-d, +d), (0.1-d, 0.1+d), **kwargs)
ax.set_ylim([y1,y2]) 
plt.gcf().subplots_adjust(bottom=0.2,left=0.2)
plt.savefig('number_of_genes_vs_constants_nice_N%i.pdf' % N,bbox_inches = "tight")


#
#import brokenaxes
#f,ax=plt.subplots()
#dummy_n_constants = np.array(n_constants,dtype=float)
#dummy_n_constants[n_constants==0] = 0.1
#
#f, bax = plt.subplots(2,1, sharex=True)
##bax = brokenaxes.brokenaxes(ylims=((np.log10(0.05),np.log10(0.2)), (np.log10(.9), np.log10(100))), hspace=.2)
##bax.plot(np.log10(n_variables),np.log10(dummy_n_constants),'ko',markersize=8,alpha=0.5)
#bax[0].plot(np.log10(n_variables),np.log10(dummy_n_constants),'ko',markersize=8,alpha=0.5)
#bax[1].plot(np.log10(n_variables),np.log10(dummy_n_constants),'ko',markersize=8,alpha=0.5)
#bax[1].set_ylim(np.log10(0.05),np.log10(0.2))
#bax[0].set_ylim(np.log10(.8), np.log10(100))
#bax[0].spines['bottom'].set_visible(False)
#bax[1].spines['top'].set_visible(False)
#bax[0].xaxis.tick_top()
#bax[0].tick_params(labeltop='off')  # don't put tick labels at the top
#bax[1].xaxis.tick_bottom()
#
#d = .03  # how big to make the diagonal lines in axes coordinates
## arguments to pass to plot, just so we don't keep repeating them
#kwargs = dict(transform=bax[0].transAxes, color='k', clip_on=False)
#bax[0].plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
#bax[0].plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
#kwargs.update(transform=bax[1].transAxes)  # switch to the bottom axes
#bax[1].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
#bax[1].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
#
#bax.set_xlabel('Number of genes')
#bax.set_ylabel('Number of constants')
#plt.gcf().subplots_adjust(bottom=0.2)
#
#plt.savefig('number_of_genes_vs_constants_nice_N%i.pdf' % N)


sorted_n_variables = np.array(sorted(zip(n_variables,n_constants),key=lambda x: x[0])).T
f,ax=plt.subplots(figsize=(10,4))
ind = np.arange(len(Fs))
ax.bar(ind,sorted_n_variables[0],width=0.7,color=[0.8,0,0],label='genes')
ax.bar(ind,sorted_n_variables[1],width=0.7,bottom=sorted_n_variables[0],color=[0.2,0.2,0.9],label='external\nparameters')
ax.set_xlim([-1,len(Fs)])
ax.set_ylabel('Number of nodes')
ax.set_xticks([])
ax.set_xlabel('Gene regulatory networks')
ax.legend(loc=2)
ax.text(len(Fs)/2.4,300,str(len(Fs))+' networks\n'+str(total_genes)+' genes\n'+str(sum(map(len,constantss)))+' external parameters',verticalalignment='center',horizontalalignment='left')
plt.savefig('summary_all_networks_N%i.pdf' % N,bbox_inches = "tight")


## proportion of increasing vs decreasing fcts
res = []
all_res = []
for i in range(len(Fs)):
    F = Fs[i]
    res.append([])
    for f in F:
        if len(f)==0: #happens if actual degree_f > max_degree
            continue
        dummy = can.is_monotonic(f,True)[1]
        res[-1].append(dummy)
        all_res.extend(dummy)
print(pd.Series(all_res).value_counts())





## Which variables occur in many networks?
variabless_simple = [[el.lower().replace('_','').replace('.','').replace('kappa','k') for el in variables] for variables in variabless]


all_variables = []
for i in range(N):
    #all_variables.extend(list(map(str.lower,variabless[i])))
    all_variables.extend(variabless_simple[i])
all_variables = np.array(all_variables)
print(pd.Series(all_variables).value_counts()[:50])
common_variables = list(pd.Series(all_variables).value_counts()[:20].index)

## In-degree and Out-degree distribution, also for common variables
all_degrees = []
for i in range(N):
    all_degrees.extend(degrees[i][:len(variabless[i])])
all_degrees = np.array(all_degrees)

all_outdegrees = []
outdegrees = []
all_outdegrees_prop = []
for i in range(N):
    outdegrees.append([0 for j in range(len(variabless[i]))])
    for regulators in Is[i][:len(variabless[i])]:
        for regulator in regulators:
            try:
                outdegrees[i][regulator] += 1
            except IndexError: #if regulator is a constant
                pass
    all_outdegrees.extend(outdegrees[i])
    mean_outdegree = np.mean(outdegrees[i])
    all_outdegrees_prop.extend([el/mean_outdegree for el in outdegrees[i]])
all_outdegrees = np.array(all_outdegrees)
all_outdegrees_prop = np.array(all_outdegrees_prop)

outdegrees_common_variables = [np.mean(all_outdegrees[all_variables==common_variables[i]]) for i in range(20)]
outdegrees_prop_common_variables = [np.mean(all_outdegrees_prop[all_variables==common_variables[i]]) for i in range(20)]

dummy = pd.Series(all_outdegrees).value_counts()
all_outdegrees_index = np.array(dummy.index)
all_outdegrees_values = np.array(dummy)
f,ax = plt.subplots()
for outdegree,count in zip(all_outdegrees_index,all_outdegrees_values):
    if outdegree>0:
        ax.loglog([outdegree],[count],'ko')
ax.set_xlabel('out-degree')
ax.set_ylabel('total count')
plt.savefig('distribution_outdegree_N%i.pdf' % len(Fs),bbox_inches = "tight")

dummy = pd.Series(all_degrees).value_counts()
all_degrees_index = np.array(dummy.index)
all_degrees_values = np.array(dummy)
f,ax = plt.subplots()
for degree,count in zip(all_degrees_index,all_degrees_values):
    if degree>0:
        ax.loglog([degree],[count],'ko')
ax.set_xlabel('in-degree')
ax.set_ylabel('total count')
plt.savefig('distribution_indegree_N%i.pdf' % len(Fs),bbox_inches = "tight")

#in degree and out degree in one plot
f,ax = plt.subplots()
ax.loglog(all_degrees_index,all_degrees_values/sum(all_degrees_values),'ro',label='in-degree')
#ax.loglog(all_outdegrees_index[all_outdegrees_index>0],all_outdegrees_values[all_outdegrees_index>0]/sum(all_outdegrees_values[all_outdegrees_index>0]),'k*',label='out-degree')
ax.loglog(all_outdegrees_index[all_outdegrees_index>0],all_outdegrees_values[all_outdegrees_index>0]/sum(all_outdegrees_values),'k*',label='out-degree')
ax.set_xlabel('degree')
ax.set_ylabel('proportion')
ax.legend(loc='best')
plt.savefig('distribution_in_and_outdegree_N%i.pdf' % len(Fs),bbox_inches = "tight")

#in degree and out degree in one plot, nice
f,ax = plt.subplots()
ax.loglog(all_degrees_index,all_degrees_values/sum(all_degrees_values),'ro',label='in-degree')
all_outdegrees_index_dummy = np.array(all_outdegrees_index,dtype=float)
all_outdegrees_index_dummy[all_outdegrees_index==0] = 0.5
#ax.loglog(all_outdegrees_index[all_outdegrees_index>0],all_outdegrees_values[all_outdegrees_index>0]/sum(all_outdegrees_values[all_outdegrees_index>0]),'k*',label='out-degree')
ax.loglog(all_outdegrees_index_dummy,all_outdegrees_values/sum(all_outdegrees_values),'k*',label='out-degree')
(x1,x2) = ax.get_xlim()
(y1,y2) = ax.get_ylim()
y1=1e-4
y2=1
ax.set_xlabel('degree')
ax.set_ylabel('proportion')
ax.legend(loc='best')
ax.set_xticks([0.5,1,2,5,10,20,50])
ax.set_xticklabels(list(map(str,[0,1,2,5,10,20,50])))
ax.spines['bottom'].set_visible(False)
ax.plot([x1,0.5],[y1,y1],'k',lw=1,clip_on=False)
ax.plot([1,x2],[y1,y1],'k',lw=1,clip_on=False)
print(x1,x2,y1,y2)
minorxticks = np.array(ax.get_xticks(minor=True))
minorxticks = minorxticks[minorxticks>1]
minorxticks = minorxticks[minorxticks<x2]
ax.set_xticks(minorxticks,minor=True)
ax.set_xlim([x1,x2]) 
ax.set_ylim([y1,y2])
#ax2 = ax.twinx()
#ax2.loglog(all_outdegrees_index_dummy,all_outdegrees_values,'k*')
#ax2.set_ylabel('total count')
#ax2.spines['bottom'].set_visible(False)
#ax2.xticks['bottom'].set_visible(False)
#ax.set_ylim([1e-4,1])
plt.savefig('distribution_in_and_outdegree_N%i_nice.pdf' % len(Fs),bbox_inches = "tight")





## aggregate all functions with k inputs
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
for i,F in enumerate(Fs):
    for f in F:
        if len(f)==0: #happens if actual degree_f > max_degree
            continue
        (NONDEGENERATED,monotonic) = can.is_monotonic(f,True)
        nrs = associate_number_monotonic(*monotonic)
        n_essential=sum([el!=-2 for el in nrs])
        n = len(nrs)
        for el in nrs:
            res.append([i,n_essential,n,el])
res = np.array(res)

color_unknown = 'black'#[0.5,0.5,0.5]
color_neg = 'orange'#[0.7,0.7,1]
color_pos = 'blue'#[1,0.7,0.7]

n_max = 7
colors = [color_pos,color_neg,color_unknown]
f,ax = plt.subplots()
x = np.arange(1,n_max+1)
res_agg = []
for i,value in enumerate([0,1,-1]):
    y = np.array([sum(res[res[:,1]==n,3]==value) for n in range(1,n_max+1)])
    res_agg.append(y)
    ax.semilogy(x[y>0],y[y>0] ,'x-',color=colors[i])
ax.legend(['increasing','decreasing','not monotonic'],loc='best')
res_agg = np.array(res_agg)

f,ax = plt.subplots()
res_agg_prop = res_agg/sum(res_agg,0)
for i,label in enumerate(['increasing','decreasing','not monotonic']):
    ax.bar(x,res_agg_prop[i,:],bottom=np.sum(res_agg_prop[:i,:],0),color=colors[i],alpha=0.3)
#ax.legend(['activation','inhibition','conditional'],bbox_to_anchor=(1.05, 0.65))
#ax.legend(['activation','inhibition','conditional'],loc=8,ncol=1)
ax.legend(['positive','negative','conditional'],loc=8,ncol=1,title='type of regulation')
ax.set_ylim([0,1])
ax.set_xticks(list(range(1,n_max+1)))
ax.set_xlabel('number of essential inputs')
ax.set_ylabel('proportion')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('activators_vs_inhibitors_vs_degree_N%i.pdf' % len(Fs),bbox_inches = "tight")


## enumerate all Boolean functions in 2 variables
res = []
n=2
numvalues = 2**n
array_numvalues = 2**np.arange(2**n-1,-1,-1)
for i,F in enumerate(Fs):
    for f in F:
        if len(f)==numvalues and can.is_degenerated(f)==False:
            res.append(np.dot(array_numvalues,f))
dict_names = {1: r'$x \wedge y$',2: r'$x \wedge \neg y$',4: r'$x \wedge \neg y$'}

dummy = pd.value_counts(res)
for i in dummy.index:
    f = can.dec2bin(i,numvalues)
    print(i,dummy[i],f,can.get_canalizing_depth_inputs_outputs_corefunction(f))
    


