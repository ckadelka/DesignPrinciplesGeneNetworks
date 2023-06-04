#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 10:19:13 2020

@author: ckadelka
"""

#Imports





import numpy as np
import itertools
import os
import canalizing_function_toolbox_v13 as can
import pickle

def find_all_indices(arr,el):
    '''
    Given a list arr, this function returns a list of all the indices i where arr[i]==el.
    If el not in arr, it raises a ValueError.
    '''
    res=[]
    for i,a in enumerate(arr):
        if a==el:
            res.append(i)
    if res==[]:
        raise ValueError('The element is not in the array at all')
    return res

def text_to_BN(folder,textfile,separator_var_func="=",original_not="NOT",original_and="AND",original_or="OR",new_not=" not ",new_and=" and ",new_or=" or ",max_degree=15,TREATMENT_OF_CONSTANTS=1,max_N=10000):
    '''
    This function takes as input a textfile in directory folder, 
    where each line describes the regulation of one gene, specified to the left, 
    separated by 'separator_var_func' from its update rules to the right.
    
    The function outputs a Boolean network model, specified by
        F: a list of N+C lists of length 2^(n_i) 
        where N is the number of genes, C the number of constants (i.e., external parameters) 
        and n_i the number of regulators per gene (n_i=1 for all constants by definition).
        F describes the update rules of the BN, stored as the right-hand side of a truth table
        
        I: a list of N lists of length n_i where N is the number of genes, n_i the number of regulators per gene. 
        I describes the regulators of each gene (as indices 0, 1, ..., n-1), i.e., the wiring diagram of the BN (the adjacency matrix of a directed graph)
        
        degree = [n_1, n_2, ..., n_N]: a list of length N, describing the number of regulators of each gene
        
        var: a list of the names of the N genes
        
        constants: a list of the names of the C constants or external parameters
    
    Inputs:
        original_not, original_and, original_or: Boolean logic operators used in the updated rules of the provided text file
        
        new_not, new_and, new_or: Boolean logic operators to be used (default: Python interpretable Boolean logic operator)
        
        max_degree: (default 15) update rules with more than max_degree regulators are not computed, instead an empty list [] is added as a placeholder
        
        max_N: (default 10,000) text files with more than max_N rows yield an AssertionError
            
        TREATMENT_OF_CONSTANTS: Ternary choice:
            0: constants (i.e., external parameters) are not added to the BN, yields a BN that cannot be dynamically evaluated and causes errors unless only the degree distribution and update rules are studied,
            1: (default) constants are added as self-regulatory nodes into the network, which is then by definition not strongly-connected,
            2: (not implemented yet) multiple models are returned, one for each combination of constants, the constants are not included as nodes but instead the update rules are simplified

    Example of an input file:
        A = NOT B
        B = A OR C
        C = E OR (A AND (NOT B))
    
    Output with TREATMENT_OF_CONSTANTS==1 (default):
        F = [[1,0],
             [0,1,1,1],
             [0,1,0,1,1,1,0,1],
             [0,1]]
        I = [[1],
             [0,2],
             [0,1,3],
             [3]]
        degree = [1,2,3,1]
        var = ['A','B','C']
        constants = ['E']
        
    Output with TREATMENT_OF_CONSTANTS==0:
        F = [[1,0],
             [0,1,1,1],
             [0,1,0,1,1,1,0,1]]
        I = [[1],
             [0,2],
             [0,1,3]]
        degree = [1,2,3]
        var = ['A','B','C']
        constants = ['E']    
    '''
    
    f = open(folder+textfile,'r')
    text = f.read()
    text = text.replace('\t',' ').replace('(',' ( ').replace(')',' ) ')
    tvec = text.splitlines()
    f.close()
    
    #Deleting empty lines
    while '' in tvec:
        tvec.remove('')
        
    n=len(tvec)
    assert n<=max_N,'n='+str(n)+' > max_N='+str(max_N)
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
        
    F = []
    for i in range(n):
        f = np.array([],dtype=int)
        if degree[i]<=max_degree:
            X = list(itertools.product([0, 1], repeat = degree[i]))
            for j in range(2**degree[i]):
                x = X[j] #x is used "implicitly" in the next line!!
                f = np.append(f,eval(tvec_mod[i])%2) #x is used here "implicitly"
        F.append(f)
        
    
    if TREATMENT_OF_CONSTANTS==1:
        for i in range(len(constants)):
            F.append(np.array([0,1]))
            I.append(np.array([len(var)+i]))
            degree.append(1)
    assert TREATMENT_OF_CONSTANTS in [0,1],'TREATMENT_OF_CONSTANTS must be 0 or 1 (default). TREATMENT_OF_CONSTANTS==2, yielding 2^C models for each combination of inputs to the C constants is not yet implemented.'
        
    return F, I, degree, var, constants

def load_tabular_model(folder,textfile,max_N=10000,TREATMENT_OF_CONSTANTS=1):
    '''
    This function takes as input a pickled file in directory folder, 
    which contains all the information about a BN, which was provided in the literature in tabular format already.
    
    The function outputs a Boolean network model, specified by
        F: a list of N+C lists of length 2^(n_i) 
        where N is the number of genes, C the number of constants (i.e., external parameters) 
        and n_i the number of regulators per gene (n_i=1 for all constants by definition).
        F describes the update rules of the BN, stored as the right-hand side of a truth table
        
        I: a list of N lists of length n_i where N is the number of genes, n_i the number of regulators per gene. 
        I describes the regulators of each gene (as indices 0, 1, ..., n-1), i.e., the wiring diagram of the BN (the adjacency matrix of a directed graph)
        
        degree = [n_1, n_2, ..., n_N]: a list of length N, describing the number of regulators of each gene
        
        var: a list of the names of the N genes
        
        constants: a list of the names of the C constants or external parameters    
    '''
    
    f = open(folder+textfile,'rb')
    [F, I, var, constants] = pickle.load(f)
    f.close()
    assert max_N>=len(var)
    degree = [len(el) for el in I]
    
    if TREATMENT_OF_CONSTANTS==1:
        for i in range(len(constants)):
            F.append([0,1])
            I.append([len(var)+i])
            degree.append(1)
    
    I = [np.array(el) for el in I]
                
    return F, I, degree, var, constants
    

def load_database(folders,separator_var_func="=",original_not="NOT",original_and="AND",original_or="OR",new_not=" not ",new_and=" and ",new_or=" or ",max_degree=15,max_N=10000):
    '''
    This function takes as input a list of directories, in which it searches for Boolean networks stored as tabular models or in text format.
    
    The function outputs a list of all M Boolean network models it was able to convert to a standardized format.
    In this format, each BN is specified by
        F: a list of N+C lists of length 2^(n_i) 
        where N is the number of genes, C the number of constants (i.e., external parameters) 
        and n_i the number of regulators per gene (n_i=1 for all constants by definition).
        F describes the update rules of the BN, stored as the right-hand side of a truth table
        
        I: a list of N lists of length n_i where N is the number of genes, n_i the number of regulators per gene. 
        I describes the regulators of each gene (as indices 0, 1, ..., n-1), i.e., the wiring diagram of the BN (the adjacency matrix of a directed graph)
        
        degree = [n_1, n_2, ..., n_N]: a list of length N, describing the number of regulators of each gene
        
        degree_essential = [e_1, e_2, ..., e_N]: a list of length N, describing the number of ESSENTIAL regulators of each gene, 
        where a regulator is essential if has an impact on the function
        
        var: a list of the names of the N genes
        
        constants: a list of the names of the C constants or external parameters
        
    The actual outputs of the function are:
        Fs: a list of M F's
        
        Is: a list of M I's
        
        degrees: a list of M degree's
        
        degrees_essential: a list of M degree_essential's
            
        variabless: a list of M var's
        
        constantss: a list of M constants's 
        
        models_loaded: a list of M filenames of the models that were successfully converted
        
        models_not_loaded: a list of filenames of the models that were NOT successfully converted
    
    Optional inputs:
        original_not, original_and, original_or: Boolean logic operators used in the updated rules of the provided text file
        
        new_not, new_and, new_or: Boolean logic operators to be used (default: Python interpretable Boolean logic operator)
        
        max_degree: (default 15) update rules with more than max_degree regulators are not computed, instead an empty list [] is added as a placeholder
        
        max_N: (default 10,000) text files with more than max_N rows yield an AssertionError
            
        TREATMENT_OF_CONSTANTS: Ternary choice:
            0: constants (i.e., external parameters) are not added to the BN, yields a BN that cannot be dynamically evaluated and causes errors unless only the degree distribution and update rules are studied,
            1: (default) constants are added as self-regulatory nodes into the network, which is then by definition not strongly-connected,
            2: multiple models are returned, one for each combination of constants, the constants are not included as nodes but instead the update rules are simplified
    '''

    Fs,Is,degrees,variabless,constantss,degrees_essential = [],[],[],[],[],[]
    models_loaded,models_not_loaded = [],[]
    
    for folder in folders:
        for fname in os.listdir(folder):
            if fname.endswith('tabular.txt'): #first check if it is a model that is already available in tabular form
                try:
                    textfile = fname
                    F,I,degree,variables, constants = load_tabular_model(folder,fname,max_N=max_N)
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
                    F,I,degree,variables, constants = text_to_BN(folder,textfile,max_degree=max_degree,max_N=max_N)
                    print(textfile,'converted')
                except:
                    models_not_loaded.append(textfile)
                    print()
                    print(textfile,'failed')
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
            degrees_essential.append([can.get_number_essential_variables(f) for f in F])
            variabless.append(variables)
            constantss.append(constants)
    return [Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded]

## Similarity of networks
def jaccard_similarity(list1, list2):
    '''
    This function computes the Jaccard similarity in [0,1] of two lists, 
    which is defined as the ratio of the size of the intersection over the size of the union.
    
    Note: The lists are turned into sets meaning duplicate entries in a list do not change the outcome.
    '''
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

## Similarity of networks
def overlap_similarity(list1, list2):
    '''
    This function computes the Jaccard similarity in [0,1] of two lists, 
    which is defined as the ratio of the size of the intersection over the size of the union.
    
    Note: The lists are turned into sets meaning duplicate entries in a list do not change the outcome.
    '''
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / min(len(list1) , len(list2))

def exclude_similar_models(Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,similarity_threshold=0.8,USE_JACCARD=True,models_to_keep=[],models_to_exclude_manually_because_similar_from_same_PID=[]):
    '''
    This function takes as input a list of Boolean network models and detects similar models 
    by comparing the Jaccard similarity of the variables of two models.
    
    Inputs: 
        The first 7 inputs come from load_database
                
        jaccard_similarity_threshold: threshold to be used to call two models similar.
        
        models_to_keep (optional): a list of models that describes which duplicate to keep (if multiple duplicates are described in this list, the first one is kept)
        
        models_to_exclude_manually_because_similar_from_same_PID (optional): a list of models that should be deleted even if not above the similarity threshold
        
    Outputs:
        The 7 inputs describing the list of Boolean network models are returned 
        but similar models have been deleted so that for each cluster of similar networks only one is kept.
        
        N: number of models are excluding similar models.
        
        models_excluded: filenames of the excluded models.
    '''
    
    N = len(Fs)
    variabless_simple = [[el.lower().replace('_','').replace('.','').replace('kappa','k').replace('-','') for el in variables] for variables in variabless]
    
    similarity_function = jaccard_similarity if USE_JACCARD else overlap_similarity
    
    sim = np.zeros((N,N))
    similar_networks = []
    dict_similar_networks = dict()
    similar_network_sets = []
    count_clusters = 0
    for i in range(N):
        for j in range(i+1,N):
            sim[i,j] = similarity_function(list(map(str.lower,variabless_simple[i])),list(map(str.lower,variabless_simple[j])))
            sim[j,i] = sim[i,j]
            if sim[i,j]>=similarity_threshold:
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
    similar_sets = [[models_loaded[i] for i in similar_network_set] for similar_network_set in similar_network_sets]
    indices_to_exclude = []
    for el,el_name in zip(similar_network_sets,similar_sets):
        dummy = []
        counter = 0
        for model_name in el_name:
            try:
                dummy.append(models_to_keep.index(model_name))
                counter+=1
            except ValueError:
                dummy.append(np.nan)
        if counter>0:
            index_to_keep = int(np.nanargmin(dummy))
            indices_to_exclude.extend(el[:index_to_keep])
            indices_to_exclude.extend(el[index_to_keep+1:])
        else:
            indices_to_exclude.extend(el[1:])
    for model_name in models_to_exclude_manually_because_similar_from_same_PID:
        try:
            index = models_loaded.index(model_name)
            indices_to_exclude.append(index)
        except ValueError:
            continue
    indices_to_exclude.sort()
    indices_to_exclude.reverse()
    
    Fs_copy = [el[:] for el in Fs]
    Is_copy = [el[:] for el in Is]
    degrees_copy = degrees[:]
    degrees_essential_copy = degrees_essential[:]
    variabless_copy = variabless[:]
    constantss_copy = constantss[:]
    models_loaded_copy = models_loaded[:]
    
    models_excluded=[]
    for index in indices_to_exclude:
        Fs_copy.pop(index)
        Is_copy.pop(index)
        degrees_copy.pop(index)
        degrees_essential_copy.pop(index)
        variabless_copy.pop(index)
        constantss_copy.pop(index)
        models_excluded.append(models_loaded[index])
        models_loaded_copy.pop(index)
        N-=1
    models_excluded.reverse()
    return (Fs_copy,Is_copy,degrees_copy,degrees_essential_copy,variabless_copy,constantss_copy,models_loaded_copy,models_excluded,similar_sets)
