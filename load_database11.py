#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 10:19:13 2020

@author: ckadelka
"""

#Imports

'''Mlc = NOT Glucose
HilE = NOT Mlc
HilD = SirA_BarA OR HilD OR HilC OR RtsA OR Fur AND (NOT CsrA) AND (NOT HilE)
HilC = HilD OR HilC OR RtsA
RtsA = HilD OR HilC OR RtsA
HilA = HilD AND (HilC OR RtsA) AND (NOT PhoP) AND ((NOT H_NS) OR (H_NS AND IHF) OR SirA)
IHF = Osmolarity AND Stationary phase
SirA_BarA = Osmolarity AND (NOT Glucose)
CsrBC = SirA_BarA
CsrA = NOT CsrBC
Fur = Iron
H_NS = NOT (PhoP AND Fur AND SlyA AND HilD AND IHF)
PhoP = NOT Magnesium
SlyA = (NOT Osmolarity) AND (NOT Calcium)
SsrAB = OmpR AND (((NOT H_NS) AND (HilD OR SlyA OR PhoP)) OR (H_NS AND (HilD OR (PhoP AND SlyA))))
Fis = Stationary_phase
EnvZ = Osmolarity
OmpR = NOT EnvZ
YfhA = NOT EnvZ
MviA = H_NS AND (NOT PhoP)
RcsB = YfhA AND (NOT MviA)
PmrA = PhoP
SciS = (RcsB AND PmrA) OR (RcsB AND PmrA AND YfhA) AND (NOT SsrAB)
VrgS = PmrA
SciG = RcsB AND PmrA
SPI_1 = HilA
SPI_2 = SsrAB AND (SlyA or Fis)
TSSS = SciS AND SciG AND VrgS'''



import numpy as np
import itertools
import os

try:
    import server.canalizing_function_toolbox_v1_9 as can
except:
    import canalizing_function_toolbox_v1_9 as can


import pickle
import traceback

def find_all_indices(array,el):
    res=[]
    for i,a in enumerate(array):
        if a==el:
            res.append(i)
    if res==[]:
        raise ValueError('The element is not in the array at all')
    return res

def text_to_BN(folder,textfile,separator_var_func="=",original_not="NOT",original_and="AND",original_or="OR",new_not=" not ",new_and=" and ",new_or=" or ",max_degree=15,TREATMENT_OF_CONSTANTS=1,max_n=10000,DONT_LOAD_IFMAX_DEGREE=False):
    '''TREATMENT_OF_CONSTANTS: Ternary choice, 
    0: constants are not added to the BN, yields a BN that cannot be dynamically evaluated and causes errors unless only the degree distribution and update rules are studied,
    1: constants are added as self-regulatory nodes into the network, which is then by definition not strongly-connected,
    2: multiple models are returned, one for each combination of constants, the constants are not included as nodes but instead the update rules are simplified'''

    f = open(folder+textfile,'r')
    text = f.read()
    f.close()
    return string_to_BN(text, separator_var_func=separator_var_func, original_not=original_not, original_and=original_and, original_or=original_or, new_not=new_not, new_and=new_and, new_or=new_or, max_degree=max_degree, TREATMENT_OF_CONSTANTS=1, max_n=max_n, DONT_LOAD_IFMAX_DEGREE=DONT_LOAD_IFMAX_DEGREE)

def string_to_BN(text,separator_var_func="=",original_not="NOT",original_and="AND",original_or="OR",new_not=" not ",new_and=" and ",new_or=" or ",max_degree=15,TREATMENT_OF_CONSTANTS=1,max_n=10000,DONT_LOAD_IFMAX_DEGREE=False):
    text = text.replace('\t',' ').replace('(',' ( ').replace(')',' ) ')
    tvec = text.splitlines()
    
    #Deleting empty lines
    while '' in tvec:
        tvec.remove('')
        
    n=len(tvec)
    assert n<=max_n,'n='+str(n)+' > max_n='+str(max_n)
    var=["" for i in range(n)]
    
    #Determining Variables, var contains the order of variables (may differ from original order)
    for i in range(n):
        var[i]=tvec[i][0:tvec[i].find(separator_var_func)].replace(' ','')

    constants_and_variables = []
    for line in tvec:
        linesplit = line.split(' ')
        for el in linesplit:
            if el not in ['(',')','+','*','1',separator_var_func,original_not,original_and,original_or,'',' ']:
                constants_and_variables.append(el)
    constants = list(set(constants_and_variables)-set(var))
        
    dict_variables_and_constants = dict({original_not:new_not,original_and:new_and,original_or:new_or})
    dict_variables_and_constants.update(dict(list(zip(var,["x[%i]" % i for i in range(len(var))]))))
    dict_variables_and_constants.update(list(zip(constants,["x[%i]" % i for i in range(len(var),len(set(constants_and_variables)))]))) #constants at end

    for i,line in enumerate(tvec):
        linesplit = line.split(' ')
        for ii,el in enumerate(linesplit):
            if el not in ['(',')','+','*','1',separator_var_func,new_not.strip(' '),new_and.strip(' '),new_or.strip(' '), '',' ']:
                linesplit[ii] = dict_variables_and_constants[el]
        tvec[i] = ' '.join(linesplit)
    #       
    for ind in range(n):
        tvec[ind]=tvec[ind][tvec[ind].find(separator_var_func)+len(separator_var_func):]
        #tvec[ind]="("+tvec[ind]+")"

    #create I, the list of essential variables
    I = []
    tvec_mod = []
    for i in range(n):
        indices_open = find_all_indices(tvec[i],'[')
        indices_end = find_all_indices(tvec[i],']')
        dummy = np.sort(np.array(list(map(int,list(set([tvec[i][(begin+1):end] for begin,end in zip(indices_open,indices_end)]))))))
        I.append( dummy )
        dict_dummy = dict(list(zip(dummy,list(range(len(dummy))))))
        tvec_dummy = tvec[i][:]
        for el in dummy:
            tvec_dummy = tvec_dummy.replace('[%i]' % el,'[%i]' % dict_dummy[el]) #needs to be done with an ascending order in dummy
        tvec_mod.append(tvec_dummy)
    
    degree = list(map(len,I))
    #assert max(degree)<=max_degree,'max_degree='+str(max(degree))+' > max_degree='+str(max_degree)+', n='+str(n)
    
    F = []
    for i in range(n):
        f = np.array([],dtype=int)
        if degree[i]<=max_degree:
            X = list(itertools.product([0, 1], repeat = degree[i]))
            for j in range(2**degree[i]):
                x = X[j]
                f = np.append(f,can.eval_expr(tvec_mod[i], x)%2)
        elif DONT_LOAD_IFMAX_DEGREE:
            raise Exception("[WARN] %d.th degree=%d is grater than max_degree=%d"%(i,degree[i],max_degree))

        F.append(f)
        
    
    if TREATMENT_OF_CONSTANTS==1:
        for i in range(len(constants)):
            F.append(np.array([0,1]))
            I.append(np.array([len(var)+i]))
            degree.append(1)
        
    return F, I, degree, var, constants


def load_tabular_model(folder,textfile,max_n=10000,TREATMENT_OF_CONSTANTS=1):
    f = open(folder+textfile,'rb')
    [F, I, var, constants] = pickle.load(f)
    f.close()
    assert max_n>=len(var)
    degree = [len(el) for el in I]
    
    if TREATMENT_OF_CONSTANTS==1:
        for i in range(len(constants)):
            F.append([0,1])
            I.append([len(var)+i])
            degree.append(1)
    
    I = [np.array(el) for el in I]
                
    return F, I, degree, var, constants
    

def load_database(folders,separator_var_func="=",original_not="NOT",original_and="AND",original_or="OR",new_not=" not ",new_and=" and ",new_or=" or ",max_degree=15,max_n=10000,DONT_LOAD_IFMAX_DEGREE=False):
    Fs,Is,degrees,variabless,constantss,degrees_essential = [],[],[],[],[],[]
    models_loaded,models_not_loaded = [],[]
    
    for folder in folders:
        for fname in os.listdir(folder):
            if fname.endswith('tabular.txt'): #first check if it is a model that is already available in tabular form
                try:
                    textfile = fname
                    F,I,degree,variables, constants = load_tabular_model(folder,fname,max_n=max_n)
                    print(textfile,'converted')
                except:
                    models_not_loaded.append(textfile)
                    print()
                    print(textfile,'failed')
                    print()
                    continue
            elif fname.endswith('.txt'):
                try:
                    textfile = fname
                    F,I,degree,variables, constants = text_to_BN(folder,textfile,max_degree=max_degree,max_n=max_n,DONT_LOAD_IFMAX_DEGREE=DONT_LOAD_IFMAX_DEGREE)
                    print(textfile,'converted')
                except:
                    models_not_loaded.append(textfile)
                    print()
                    print(textfile,'failed')
                    #traceback.print_exc()
                    print()
                    continue
            else:
                continue
            #if len(constants)>0:
            #    print(textfile,'has %i constants' % len(constants))
            
            models_loaded.append(textfile)            
            Fs.append(F)
            Is.append(I)
            degrees.append(degree)
            degrees_essential.append([can.nr_essential_variables(f) for f in F])
            variabless.append(variables)
            constantss.append(constants)
    return [Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded]

## Similarity of networks
def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

def exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,jaccard_similarity_threshold=0.8):
    variabless_simple = [[el.lower().replace('_','').replace('.','').replace('kappa','k') for el in variables] for variables in variabless]
    
    sim = np.zeros((N,N))
    similar_networks = []
    dict_similar_networks = dict()
    similar_network_sets = []
    count_clusters = 0
    for i in range(N):
        for j in range(i+1,N):
            #sim[i,j] = jaccard_similarity(list(map(str.lower,variabless[i])),list(map(str.lower,variabless[j])))
            sim[i,j] = jaccard_similarity(list(map(str.lower,variabless_simple[i])),list(map(str.lower,variabless_simple[j])))
            sim[j,i] = sim[i,j]
            if sim[i,j]>jaccard_similarity_threshold:
                similar_networks.append([i,j,sim[i,j],models_loaded[i],models_loaded[j],len(Fs[i]),len(Fs[j])])
                try:
                    cluster_id = dict_similar_networks[i]
                except KeyError:
                    try:
                        cluster_id = dict_similar_networks[j]
                    except KeyError:
                        cluster_id = count_clusters
                        count_clusters += 1
                        similar_network_sets.append(set())
                dict_similar_networks.update({i:cluster_id})
                dict_similar_networks.update({j:cluster_id})
                similar_network_sets[cluster_id].add(i)
                similar_network_sets[cluster_id].add(j)
        
    similar_network_sets = list(map(list,similar_network_sets))
    indices_to_exclude = []
    for el in similar_network_sets:
        indices_to_exclude.extend(el[1:])
    indices_to_exclude.sort()
    indices_to_exclude.reverse()
    
    models_excluded=[]
    for index in indices_to_exclude:
        Fs.pop(index)
        Is.pop(index)
        degrees.pop(index)
        degrees_essential.pop(index)
        variabless.pop(index)
        constantss.pop(index)
        models_loaded.pop(index)
        N-=1
        models_excluded.append(models_loaded[index])
    return (Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,N,models_excluded)
