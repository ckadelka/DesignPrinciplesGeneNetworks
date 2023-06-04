#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:01:32 2019

@author: ckadelka
"""

##Imports

import numpy as np
import matplotlib.pyplot as plt
import canalizing_function_toolbox_v1_9 as can
import scipy.stats as stats
import load_database11 as db

## load the database, choose low max_n for quick results and to only look at small models
folders = ['update_rules_cell_collective/', 'update_rules_models_in_literature_we_randomly_come_across/']
max_degree = 16
max_n=1000
[Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded] = db.load_database(folders,max_degree=max_degree,max_n=max_n)
N = len(models_loaded)
jaccard_similarity_threshold = 0.8
Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded = db.exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=jaccard_similarity_threshold)

all_max_degrees = [max(degree) for degree in degrees]

## Derrida plot - some basic plots
n_variables = [len(el) for el in variabless]
n_cs = [len(el) for el in constantss]
derridas = []
for i in range(len(Fs)):
    if all_max_degrees[i]>max_degree:
        continue
    F = Fs[i]
    I = Is[i]
    len_F = len(F)
    derrida = can.derrida_value(F,I,len_F,1,nsim=1000)
    derridas.append(derrida)

n_variables_restr = [n_variables[i] for i in range(N) if all_max_degrees[i]<=max_degree]

plt.hist(derridas)
plt.scatter(n_variables_restr,derridas)
plt.scatter(n_variables_restr,derridas)
stats.spearmanr(n_variables_restr,derridas)

f,ax = plt.subplots()
ax.semilogx(n_variables_restr,derridas,'ko',alpha=0.5)
ax.set_xlabel('Number of genes')
ax.set_ylabel('Derrida value for a\n single perturbation')
(x1,x2) = ax.get_xlim()
(y1,y2) = ax.get_ylim()
r,p = stats.spearmanr(n_variables_restr, derridas)
ax.text(np.exp(np.log(x1)+0.73*(np.log(x2)-np.log(x1))),y1+0.87*(y2-y1),(r'r$=%s$' % str(np.round(r,2))) +'\n'+('p$=%s$' % str(np.round(p,2))),color='red',va='center',ha='left')
plt.gcf().subplots_adjust(bottom=0.2)
plt.savefig('size_vs_derrida_N%i.pdf' % len(Fs),bbox_inches = "tight")   


import seaborn as sns
f = sns.jointplot(np.log10(n_variables_restr),derridas)
x = np.log10(n_variables_restr)
fontsize=16
g = sns.jointplot(x=x, y=derridas, kind="reg",color="k",ratio=2,marginal_kws=dict(bins=20,color='b'))
g.ax_joint.set_xticks([1,2,3],minor=False)
minorticks = np.r_[np.log10(1*np.arange(2,10)),np.log10(10*np.arange(2,10)),np.log10(100*np.arange(2,10))]
g.ax_joint.set_xticks(minorticks,minor=True)
print(g.ax_joint.get_yticks())
g.ax_joint.set_xticklabels(['10','100','1000'],fontsize=fontsize)
g.ax_joint.set_yticklabels([str(round(el,1)) for el in g.ax_joint.get_yticks()],fontsize=fontsize)
g.ax_joint.set_xlim(np.log10(2),np.log10(400))
g.ax_joint.set_xlabel('Number of genes',fontsize=fontsize)
g.ax_joint.set_ylabel('Average size of a single\nperturbation after one update',fontsize=fontsize)
g.savefig('size_vs_derrida_N%i_nice.pdf' % len(Fs),bbox_inches = "tight")

f, ax = plt.subplots(figsize=(5, 5))
ax.set(xscale="log")
g = sns.jointplot(np.log10(n_variables_restr),derridas)
g.ax_joint.set_xscale('log')
g.ax_joint.set_xticks([1,2])
g.ax_joint.set_xlim([0.3,2.5])


sns.regplot("x", "y", data, ax=ax, scatter_kws={"s": 100})
