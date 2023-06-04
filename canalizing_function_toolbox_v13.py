#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 15:53:40 2019

@author: ckadelka

#Version 13

This toolbox allows to 

1) determine if a Boolean function is 
a) constant
b) degenerated
c) canalizing
d) k-canalizing

2) determine the 
a) canalizing depth of a Boolean function 
b) the layer structure of a canaliizing Boolean function

3) generate uniformly at random
a) non-degenerated Boolean functions
a) non-canalizing Boolean functions
c) non-canalizing non-degenerated Boolean functions
d) k-canalizing Boolean functions 
e) k-canalizing Boolean functions with a given layer structure
f) Boolean functions with exact canalizing depth k
g) Boolean functions with exact canalizing depth k with a given layer structure

4) generate uniformly at random Boolean networks with specific characterists (in-degree, canalization, strong connectedness)

5) obtain some basic estimates of the magnitude of various subclasses of Boolean functions

6) determine the 
a) absolute bias of a Boolean function
b) average sensitivity of a Boolean function
"""

#13:  proper documentation added, deleted functions that became obsolete
#1.9: new functionality added: calculate feed forward loops and feedback loops
#1.5: added is_kset_canalizing
#1.4: fixed issue with k==0 and EXACT_DEPTH==True in random_BN
#1.3: Python3.7 compatible, kis passed to random_BN can also be [0] for random networks
#1.2: added functionality to randomly create and analyze Boolean networks based on random or canalizing functions
#1.1: fixed a couple issues in is_k_canalizing, is_k_canalizing_return_inputs_outputs_corefunction and get_layerstructure_given_canalizing_outputs_and_corefunction

##Imports

import numpy as np
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import random

import sympy
import pandas as pd

from collections import defaultdict
from matplotlib import colors as mcolors

try:
    import cana.boolean_node
    LOADED_CANA=True
except:
    LOADED_CANA=False

## 0) Basics, helper functions    
def tobin(x):
    '''returns the binary representation (in array form) of a decimal number'''
    return tobin(x//2) + [x%2] if x > 1 else [x]

def dec2bin(x,n=[]):
    '''returns the binary representation (in array form) of a decimal number.
    Input can be an array itself, in which case each decimal number in the 
    array is separately converted into a binary array.
    The second input n, if specified, describes the number of positions in the
    binary representation array for each number.
    Example: dec2bin(10)=dec2bin(10,4)=[1,0,1,0]
             dec2bin(10,6)=[0,0,1,0,1,0]'''
    if type(x) in [list,np.ndarray]:
        return [dec2bin(el,n) for el in x]
    if n==[]:
        return tobin(x)
    else:
        help=tobin(x)
        res=[0]*(n-len(help))
        res.extend(help)
        return res

def bin2dec(state):
    '''returns the decimal number representation of a binary state'''
    n = len(state)
    b = [2**i for i in range(n)]
    return sum([state[n-i-1]*b[i] for i in range(n)])

def find_all_indices(array,value):
    '''returns a list of all indices in array that equal a given value'''
    res=[]
    for i,a in enumerate(array):
        if a==value:
            res.append(i)
    if res==[]:
        raise ValueError('The value is not in the array at all.')
    return res

def edgelist_to_I(edgelist):
    '''input: an m x 2 array describing all edges (i.e., regulations), 
    with the first column describing the regulator, and the second the regulated node.
    
    outputs: 
        1. I, a list of lists where 
        2. var, a list of all variables that show up as regulators and/or regulated nodes.
    '''    
    regulators = np.array(edgelist)[:,0]
    targets = np.array(edgelist)[:,1]
    var = list(set(regulators)|set(targets))
    n_var = len(var)
    dict_var = dict(zip(var,range(n_var)))
    I = [[] for _ in range(n_var)]
    for i in range(len(regulators)):
        I[dict_var[targets[i]]].append(dict_var[regulators[i]])
    return I,var

def bool_to_poly(f,left_side_of_truth_table=[]):
    '''
    This function transforms a Boolean function from truth table format to polynomial format.
    The polynomial is in non-reduced disjunctive normal form (DNF).

    inputs:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
        left_side_of_truth_table (optional for speed-up): the left-hand side of the Boolean truth table of size 2^n x n
    
    output: a string of the Boolean function in disjunctive normal form (DNF)
    '''
    len_f = len(f)
    n=int(np.log2(len_f))
    if left_side_of_truth_table==[]: #to reduce run time, this should be calculated once and then passed as argument
        left_side_of_truth_table = list(itertools.product([0, 1], repeat = n))
    num_values = 2**n   
    text = []
    for i in range(num_values):
        if f[i]==True:
            monomial = '*'.join([('x%i' % (j+1)) if entry==1 else ('(1-x%i)' % (j+1)) for j,entry in enumerate(left_side_of_truth_table[i])])
            text.append(monomial)
    if text!=[]:
        return ' + '.join(text)
    else:
        return '0'
            
## 1) Methods to analyze Boolean functions

def is_degenerated(f):
    '''
    This function determines if a Boolean function contains some non-essential variables.
    
    input: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
    output: TRUE if f contains non-essential variables, FALSE if all variables are essential.
    '''    
    len_f = len(f)
    n=int(np.log2(len_f))
    for i in range(n):
        dummy_add=(2**(n-1-i))
        dummy=np.arange(2**n)%(2**(n-i))//dummy_add
        depends_on_i=False
        for j in range(2**n):
            if dummy[j]==1:
                continue
            else:
                if f[j]!=f[j+dummy_add]:
                    depends_on_i=True
                    break
        if depends_on_i==False:
            return True
    return False

def get_essential_variables(f):
    '''
    This function determines the essential variables of a Boolean function 
    by testing exhaustively whether a given variable is essential.
    
    input: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
    output: a list of all indices of variables, which are essential 
    '''    
    if len(f)==0:
        return []
    len_f = len(f)
    n=int(np.log2(len_f))
    essential_variables  = list(range(n))
    for i in range(n):
        dummy_add=(2**(n-1-i))
        dummy=np.arange(2**n)%(2**(n-i))//dummy_add
        depends_on_i=False
        for j in range(2**n):
            if dummy[j]==1:
                continue
            else:
                if f[j]!=f[j+dummy_add]:
                    depends_on_i=True
                    break
        if depends_on_i==False:
            essential_variables.remove(i)
    return essential_variables 

def get_number_essential_variables(f):
    '''
    This function determines the number of essential variables of a Boolean function.
    
    input: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
    output: an integer, the number of essential variables of f  
    '''        
    return len(get_essential_variables(f)) 

def is_constant(f):
    '''
    This function checkes whether a Boolean function is constant.

    input: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
    output: TRUE if f is constant, FALSE otherwise.    
    '''
    return sum(f) in [0,len(f)]

def get_symmetry_groups(f,left_side_of_truth_table=[]):
    '''
    This function determines all symmetry groups of input variables for a Boolean function.
    Two variables x,y are in the same symmetry group if f(x,y,z1,...,z_m) = f(y,x,z1,...,z_m) for all possible inputs to the other variables z1,...,z_m

    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
        left_side_of_truth_table (optional for speed-up): the left-hand side of the Boolean truth table of size 2^n x n
    
    output: a list of m lists where m is the number of symmetry groups of f 
    and each inner list contains the indices of all variables in the same symmetry group.
    '''
    len_f = len(f)
    n=int(np.log2(len_f))
    if left_side_of_truth_table==[] or left_side_of_truth_table.shape[0]!=len_f:
        left_side_of_truth_table = np.array(list(itertools.product([0, 1], repeat=n)))
    symmetry_groups = []
    left_to_check = np.ones(n)
    for i in range(n):
        if left_to_check[i]==0:
            continue
        else:
            symmetry_groups.append([i])
            left_to_check[i]=0
        for j in range(i+1,n):
            diff = sum(2**np.arange(n-i-2,n-j-2,-1))
            for ii,x in enumerate(left_side_of_truth_table):
                if x[i]!=x[j] and x[i]==0 and f[ii]!=f[ii+diff]:
                    break
            else:
                left_to_check[j] = 0
                symmetry_groups[-1].append(j)
    return symmetry_groups

def is_canalizing(f,n=-1):
    '''
    This function determines if a Boolean function is canalizing.
    A Boolean function f(x_1,...,x_n) is canalizing if it is canalizing in at least one variable.
    A Boolean function f(x_1,...,x_n) is canalizing in x_i if f(x_1,...,x_i=a,...,x_n) = b for some a,b in [0,1] and for all x_1,...,x_{i-1},x_{i+1},...,x_n in [0,1].

    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
        n (optional for minor speed-up): number of variables of f
    
    output: TRUE if f is canalizing, FALSE otherwise.
    '''
    if type(f) == list:
        f = np.array(f)
    if n==-1:
        n=int(np.log2(len(f)))
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    Atimesf = np.dot(A,f)
    if np.any(Atimesf==desired_value):
        return True
    elif np.any(Atimesf==0):
        return True
    else:
        return False

def is_kset_canalizing(f,k,n=-1):
    '''
    This function determines if a Boolean function is k-set canalizing.
    A Boolean function f(x_1,...,x_n) is k-set canalizing 
    if there exists a set of k variables such that if this set of variables takes on certain inputs, 
    the output of f is determined, irrespective of the input to the n-k other variables.

    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
        k, an integer (meaningful integers k in {0,1,...,n})
    
        n (optional for minor speed-up): number of variables of f
    
    output: TRUE if f is k-set canalizing, FALSE otherwise.
    
    references:
        Kadelka, C., Keilty, B., & Laubenbacher, R. (2023). Collectively canalizing Boolean functions. Advances in Applied Mathematics, 145, 102475.
    '''
    if type(f) == list:
        F = np.array(f)
    if k==0:
        return is_constant(F)
    if n==-1:
        n=int(np.log2(len(f)))
    desired_value = 2**(n-k)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    Ak = []
    for i in range(2*n):
        for j in range(i+1,2*n):
            if j-i == n:
                continue
            else:
                Ak.append( np.bitwise_and(A[i,:],A[j,:]) )
                
    Ak = []                
    for indices in itertools.combinations(range(2*n),k):
        dummy = np.sum(A[np.array(indices),:],0)==k
        if sum(dummy)==desired_value:
            Ak.append(dummy)
    Ak = np.array(Ak)
    
    AktimesF = np.dot(Ak,F)
    is_there_canalization = 0 in AktimesF or desired_value in AktimesF
    return is_there_canalization
    
def get_proportion_of_collectively_canalizing_input_sets(f,k,n=-1,left_side_of_truth_table=[],verbose=False):
    '''
    A Boolean function f(x_1,...,x_n) is k-set canalizing 
    if there exists a set of k variables such that if this set of variables takes on certain inputs, 
    the output of f is determined, irrespective of the input to the n-k other variables.
    For a given k, this function computes the probability that a k-set canalizes f.

    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
        k, an integer (meaningful integers k in {0,1,...,n})
    
        n (optional for minor speed-up): number of variables of f
        
        left_side_of_truth_table (optional for speed-up): the left-hand side of the Boolean truth table of size 2^n x n

        verbose (optional bool): TRUE to print all canalizing k-sets
    
    output: TRUE if f is k-set canalizing, FALSE otherwise.
    
    references:
        Kadelka, C., Keilty, B., & Laubenbacher, R. (2023). Collectively canalizing Boolean functions. Advances in Applied Mathematics, 145, 102475.
    '''
    if type(f) == list:
        f = np.array(f)
    if k==0:
        return float(is_constant(f))
    if n==-1:
        n=int(np.log2(len(f)))
    desired_value = 2**(n-k)
    if left_side_of_truth_table == []:
        T = np.array(list(itertools.product([0, 1], repeat=n))).T
    else:
        T = np.array(left_side_of_truth_table).T  
    Tk = list(itertools.product([0, 1], repeat=k))
    A = np.r_[T,1-T]
    Ak = []                
    for indices in itertools.combinations(range(n),k):
        for canalizing_inputs in Tk:
            indices_values = np.array(indices) + n*np.array(canalizing_inputs)
            dummy = np.sum(A[indices_values,:],0)==k
            if sum(dummy)==desired_value:
                Ak.append(dummy)
                if verbose and np.dot(dummy,f) in [0,desired_value]:
                    print(indices,canalizing_inputs,indices_values,np.dot(dummy,f))
            elif verbose:
                print(indices,canalizing_inputs,sum(dummy),'a')
    Ak = np.array(Ak)
    is_there_canalization = np.in1d(np.dot(Ak,f),[0,desired_value])    
    return sum(is_there_canalization)/len(is_there_canalization)

def binom(n,k):
    import scipy.special
    return scipy.special.binom(n,k)

def get_canalizing_strength(f,left_side_of_truth_table=[]):
    '''
    This function computes the canalizing strength of a Boolean function by exhaustive enumeration (slow for functions with many variables).
    The canalizing strength is a weighted average of the 1- to (n-1)-set canalizing proportions.
    It is 0 for the least canalizing functions, Boolean parity functions (e.g., f= (x1 + x2 + ... + xn) % 2 == 0)
    and is 1 for the most canalizing non-constant functions, nested canalizing functions with one layer (e.g. f = x1 & x2 & ... & xn)

    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
        
        left_side_of_truth_table (optional for speed-up): the left-hand side of the Boolean truth table of size 2^n x n
    
    output: a float, describing the canalizing strength of f
    
    references:
        Kadelka, C., Keilty, B., & Laubenbacher, R. (2023). Collectively canalizing Boolean functions. Advances in Applied Mathematics, 145, 102475.
    '''
    nfloat = np.log2(len(f))
    n = int(nfloat)
    assert abs(n-nfloat)<1e-10, "F needs to be of length 2^n for some n>1"
    assert n>1, "Canalizing strength is only defined for Boolean functions with n>1 inputs"
    res = []
    for k in range(1,n):
        res.append(get_proportion_of_collectively_canalizing_input_sets(f,k,n,left_side_of_truth_table=left_side_of_truth_table))
    return np.mean(np.multiply(res,2**np.arange(1,n)/(2**np.arange(1,n)-1))),res

def compute_exact_kset_canalizing_proportion_for_ncf_with_specific_layerstructure(k,layerstructure_NCF):
    '''
    This function implements Theorem 3.3 in [1]. 
    It computes the exact k-set canalizing proportion for a Boolean NCF with known layer structure.
    
    input:
        k, an integer (meaningful integers k in {0,1,...,n}) where n is the number of variables of the NCF
        
        layerstructure_NCF: [k_1,..., k_r] a list of integers describing the number of variables in each layer of an NCF,
        k_i >= 1, and k_r >= 2 unless r = n = 1. 
    
    output: a float, describing the k-set canalizing proportion for the NCF with the provided layer structure
    
    references:
        [1] Kadelka, C., Keilty, B., & Laubenbacher, R. (2023). Collectively canalizing Boolean functions. Advances in Applied Mathematics, 145, 102475.
    '''
    r = len(layerstructure_NCF)
    n = sum(layerstructure_NCF)
    assert min(layerstructure_NCF) >= 1 and (layerstructure_NCF[-1]>= 2 or n==1), "each layer must contain at least one variable (the last layer at least two unless n==1)"
    magnitudes = []
    for t in range(r):
        number_of_input_sets = 0
        for c in range(1,min( k-sum(layerstructure_NCF[:t][::-2]) , layerstructure_NCF[t] )+1):
            for d in range(0,min(k-sum(layerstructure_NCF[:t][::-2])-c , sum(layerstructure_NCF[:max(0,t-1)][::-2]))+1):
                binom1 = binom( layerstructure_NCF[t] , c )
                binom2 = binom( sum(layerstructure_NCF[:max(0,t-1)][::-2]) , d )
                binom3 = binom( n-sum(layerstructure_NCF[:t+1]) , k - sum(layerstructure_NCF[:t][::-2]) - c - d)
                number_of_inputs_that_canalize_for_selected_variable_set = sum([2**( k - sum(layerstructure_NCF[:t][::-2]) - j - d ) for j in range(1,c+1)])
                number_of_input_sets += binom1 * binom2 * binom3 * number_of_inputs_that_canalize_for_selected_variable_set
        magnitudes.append(number_of_input_sets)
    #for the case where the non-canalizing output value can be reached in the evaluation process, add:
    if k >= sum(layerstructure_NCF[-1::-2]):
        magnitudes.append( binom( n-sum(layerstructure_NCF[-1::-2]), k-sum(layerstructure_NCF[-1::-2]) ) ) 
    else:
        magnitudes.append( 0 ) 
    return sum(magnitudes)/(2**k * binom(n,k))#, magnitudes


#test using this code:
# k=3; kis=[2,1,2]; 
# print(compute_kset_canalizing_proportion_for_ncf(k,kis)); 
# print(can.get_canalizing_strength(can.random_k_canalizing_with_specific_layerstructure(sum(kis),kis))[1][k-1])

def is_k_canalizing(f,k,n=-1):
    '''
    This function determines if a Boolean function is k-canalizing.
    A Boolean function f(x_1,...,x_n) is k-canalizing if it has at least k conditionally canalizing variables.
    In other words, if 
    1. f is canalizing, and if
    2. the subfunction of f when the canalizing variable takes on its non-canalizing function is canalizing, and if
    3. the subfunction of the subfunction when its canalizing variable takes on its non-canalizing function is canalizing, etc. 
    The number of such variables is the canalizing depth of a Boolean function 
    and a function with canalizing depth >= k is k-canalizing.
    
    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
        k, an integer (meaningful integers k in {0,1,...,n}).
        Note: any function is 0-canalizing. Only NCFs are n-canalizing.
    
        n (optional for minor speed-up): number of variables of f
    
    output: TRUE if f is k-canalizing, FALSE otherwise.
    
    references:
        He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth. Physica D: Nonlinear Phenomena, 314, 1-8.

        Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions: Algorithms and applications. Automatica, 146, 110630.
    '''
    if k>n:
        return False
    if k==0:
        return True
    if n==-1:
        n=int(np.log2(len(f)))
    w = sum(f) #Hamming weight of f
    if w == 0 or w == 2**n: #constant f
        return False
    if type(f) == list:
        f = np.array(f)
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    try: #is 1 canalized output for one of the variables
        index = list(np.dot(A,f)).index(desired_value)
        new_f = f[np.where(A[index]==0)[0]]
        return is_k_canalizing(new_f,k-1,n-1)
    except ValueError:
        try: #is 0 canalized output for one of the variables
            index = list(np.dot(A,1-f)).index(desired_value)
            new_f = f[np.where(A[index]==0)[0]]
            return is_k_canalizing(new_f,k-1,n-1)
        except ValueError:
            return False
    

def is_k_canalizing_return_inputs_outputs_corefunction(f,k,n,can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int)):
    '''
    This function determines if a Boolean function is k-canalizing.
    A Boolean function f(x_1,...,x_n) is k-canalizing if it has at least k conditionally canalizing variables.
    In other words, if 
    1. f is canalizing, and if
    2. the subfunction of f when the canalizing variable takes on its non-canalizing function is canalizing, and if
    3. the subfunction of the subfunction when its canalizing variable takes on its non-canalizing function is canalizing, etc. 
    The number of such variables is the canalizing depth of a Boolean function 
    and a function with canalizing depth >= k is k-canalizing.
    
    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
        k, an integer (meaningful integers k in {0,1,...,n}).
        Note: any function is 0-canalizing. Only NCFs are n-canalizing.
    
        n (optional for minor speed-up): number of variables of f
    
    output: 
        bool: TRUE if f is k-canalizing, FALSE otherwise.
        
        can_inputs: list of the first k canalizing input values of f if f is k-canalizing, otherwise list of all canalizing input values

        can_outputs: list of the first k canalized output values of f if f is k-canalizing, otherwise list of all canalized output values

        core_function: Boolean core function in n-k variables if f is k-canalizing, otherwise Boolean core function in all m>n-k non-conditionally canalizing variables
    
    references:
        He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth. Physica D: Nonlinear Phenomena, 314, 1-8.

        Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions: Algorithms and applications. Automatica, 146, 110630.
    '''
    if k==0:
        return (True,can_inputs,can_outputs,f)
    w = sum(f)
    if w == 0 or w == 2**n: #constant f
        return (False,can_inputs,can_outputs,f)
    if type(f) == list:
        f = np.array(f)
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    try: #is 1 canalized output for one of the variables
        index = list(np.dot(A,f)).index(desired_value)
        new_f = f[np.where(A[index]==0)[0]]
        return is_k_canalizing_return_inputs_outputs_corefunction(new_f,k-1,n-1,np.append(can_inputs,int(index<n)),np.append(can_outputs,1))
    except ValueError:
        try: #is 0 canalized output for one of the variables
            index = list(np.dot(A,1-f)).index(desired_value)
            new_f = f[np.where(A[index]==0)[0]]
            return is_k_canalizing_return_inputs_outputs_corefunction(new_f,k-1,n-1,np.append(can_inputs,int(index<n)),np.append(can_outputs,0))
        except ValueError:
            return (False,can_inputs,can_outputs,f)     

def is_k_canalizing_return_inputs_outputs_corefunction_order(F,k,n,can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int),can_order=np.array([],dtype=int),variables=[]):
    '''
    This function determines if a Boolean function is k-canalizing.
    A Boolean function f(x_1,...,x_n) is k-canalizing if it has at least k conditionally canalizing variables.
    In other words, if 
    1. f is canalizing, and if
    2. the subfunction of f when the canalizing variable takes on its non-canalizing function is canalizing, and if
    3. the subfunction of the subfunction when its canalizing variable takes on its non-canalizing function is canalizing, etc. 
    The number of such conditionally canalizing variables is the canalizing depth of f, 
    and a function with canalizing depth >= k is k-canalizing.
    
    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
        k, an integer (meaningful integers k in {0,1,...,n}).
        Note: any function is 0-canalizing. Only NCFs are n-canalizing.
    
        n (optional for minor speed-up): number of variables of f
    
    output: 
        bool: TRUE if f is k-canalizing, FALSE otherwise.
        
        can_inputs: list of the first k canalizing input values of f if f is k-canalizing, otherwise list of all canalizing input values

        can_outputs: list of the first k canalized output values of f if f is k-canalizing, otherwise list of all canalized output values

        core_function: Boolean core function in n-k variables if f is k-canalizing, otherwise Boolean core function in all m>n-k non-conditionally canalizing variables
    
        can_order: list of the indices of the first k conditionally canalizing variables if f is k-canalizing, 
        otherwise list of the indices of all conditionally canalizing variables  
    
    references:
        He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth. Physica D: Nonlinear Phenomena, 314, 1-8.

        Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions: Algorithms and applications. Automatica, 146, 110630.
    '''
    if k==0:
        return (True,can_inputs,can_outputs,F,can_order)
    w = sum(F)
    if w == 0 or w == 2**n: #constant F
        return (False,can_inputs,can_outputs,F,can_order)
    if type(variables)==np.ndarray:
        variables = list(variables)
    if variables == []:
        variables = list(range(n))
    if type(F) == list:
        F = np.array(F)
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    try: #is 0 canalized output for one of the variables
        index = list(np.dot(A,1-F)).index(desired_value)
        newF = F[np.where(A[index]==0)[0]]
        variable = variables.pop(index%n)
        return is_k_canalizing_return_inputs_outputs_corefunction_order(newF,k-1,n-1,np.append(can_inputs,int(index<n)),np.append(can_outputs,0),np.append(can_order,variable),variables)
    except ValueError:
        try: #is 1 canalized output for one of the variables
            index = list(np.dot(A,F)).index(desired_value)
            newF = F[np.where(A[index]==0)[0]]
            variable = variables.pop(index%n)
            return is_k_canalizing_return_inputs_outputs_corefunction_order(newF,k-1,n-1,np.append(can_inputs,int(index<n)),np.append(can_outputs,1),np.append(can_order,variable),variables)
        except ValueError:
            return (False,can_inputs,can_outputs,F,can_order)    

def find_layers(f,can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int),can_order=np.array([],dtype=int),variables=[],depth=0,number_layers=0):
    '''
    find_layers from https://github.com/ckadelka/BooleanCanalization/blob/main/find_layers.py
    
    By [1], any non-zero Boolean function has a unqiue standard monomial form, in which
    all conditionally canalizing variables are distributed into layers of importance.
    This function determines the canalizing layer format of a Boolean function.
    It is Algorithm 2 from [2]. For a fast implementation of Algorithm 2, see the original github repo.
    
    input:
        f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)

        all other optional arguments must not be specified, for recursion only
    
    output: 
        depth: integer k>=0, number of conditionally canalizing variables
            
        number_layers: integer >=0, number of different layers
        
        can_inputs: list of all k canalizing input values of f 

        can_outputs: list of all k canalized output values of f

        core_polynomial: Boolean core polynomial in all n-k non-conditionally canalizing variables 
    
        can_order: list of the indices of all conditionally canalizing variables
    
    references:
        [1] He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth. Physica D: Nonlinear Phenomena, 314, 1-8.

        [2] Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions: Algorithms and applications. Automatica, 146, 110630.
    '''
    n = int(np.log2(len(f)))
    w = sum(f)
    if w == 0 or w == 2**n: #constant f
        return (depth,number_layers,can_inputs,can_outputs,f,can_order)
    if type(variables)==np.ndarray:
        variables = list(variables)
    if variables == []:
        variables = list(range(n))
    if type(f) == list:
        f = np.array(f)
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
        
    indices1 = np.where(np.dot(A,f)==desired_value)[0]
    indices0 = np.where(np.dot(A,1-f)==desired_value)[0]
    if len(indices1)>0:
        sorted_order = sorted(range(len(indices1)),key=lambda x: (indices1%n)[x])
        inputs = (1-indices1//n)[np.array(sorted_order)]
        outputs = np.ones(len(indices1),dtype=int)
        new_canalizing_variables = []
        for index in np.sort(indices1%n)[::-1]:
            new_canalizing_variables.append(variables.pop(index))
        new_canalizing_variables.reverse()
        new_f = f[np.sort(list(set.intersection(*[] + [set(np.where(A[index]==0)[0]) for index,INPUT in zip(indices1,inputs)])))]
        return find_layers(new_f,np.append(can_inputs,inputs),np.append(can_outputs,outputs),np.append(can_order,new_canalizing_variables),variables,depth+len(new_canalizing_variables),number_layers+1)
    elif len(indices0): 
        sorted_order = sorted(range(len(indices0)),key=lambda x: (indices0%n)[x])
        inputs = (1-indices0//n)[np.array(sorted_order)]
        outputs = np.zeros(len(indices0),dtype=int)
        new_canalizing_variables = []#variables[indices0%n]
        for index in np.sort(indices0%n)[::-1]:
            new_canalizing_variables.append(variables.pop(index))
        new_canalizing_variables.reverse()
        new_f = f[np.sort(list(set.intersection(*[] + [set(np.where(A[index]==0)[0]) for index,INPUT in zip(indices0,inputs)])))]
        return find_layers(new_f,np.append(can_inputs,inputs),np.append(can_outputs,outputs),np.append(can_order,new_canalizing_variables),variables,depth+len(new_canalizing_variables),number_layers+1)
    else:
        return (depth,number_layers,can_inputs,can_outputs,f,can_order)

## 2) Put everything together to obtain canalizing depth, layer structure, canalized outputs, canalizing inputs as well as core function (could also calculate canalizing variables in future versions but I don't see a need)

if LOADED_CANA:
    def get_input_redundancy(f,n=-1):
        '''
        This function computes the input redundancy of a Boolean function, defined as in [1].
        Constant functions have an input redundancy of 1
        because none of the inputs are needed to know the output of the function.
        Parity functions (e.g., f= (x1 + x2 + ... + xn) % 2 == 0) have an input redundancy of 0
        because all inputs are always needed to know the output of the function.

        input:
            f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
        
        output: a float in [0,1] describing the normalized input redundancy of f

        references:
            [1] Marques-Pita, M., & Rocha, L. M. (2013). Canalization and control in automata networks: body segmentation in Drosophila melanogaster. PloS one, 8(3), e55946.
            
            [2] Correia, R. B., Gates, A. J., Wang, X., & Rocha, L. M. (2018). CANA: a python package for quantifying control and canalization in Boolean networks. Frontiers in physiology, 9, 1046.
        '''
        if n==-1:
            n = int(np.log2(len(f)))
        return cana.boolean_node.BooleanNode(k=n,outputs=f).input_redundancy()

    def get_edge_effectiveness(f,n=-1):
        '''
        This function computes the edge effectiveness for each regulator of a Boolean function, defined as in [1].
        Non-essential inputs have an effectiveness of 0.
        Inputs, which when flipped always flip the output of a function, have an effectiveness of 1.

        input:
            f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)

        output: a list of N floats in [0,1] describing the edge effectiveness of each regulator of f

        references:
            [1] Marques-Pita, M., & Rocha, L. M. (2013). Canalization and control in automata networks: body segmentation in Drosophila melanogaster. PloS one, 8(3), e55946.
            
            [2] Correia, R. B., Gates, A. J., Wang, X., & Rocha, L. M. (2018). CANA: a python package for quantifying control and canalization in Boolean networks. Frontiers in physiology, 9, 1046.
        '''
        if n==-1:
            n = int(np.log2(len(f)))
        return cana.boolean_node.BooleanNode(k=n,outputs=f).edge_effectiveness()

def get_canalizing_depth_inputs_outputs_corefunction(f):
    '''obsolete, included for backward compatability, use find_layers(f)'''
    n = int(np.log2(len(f)))
    (NESTED_CANALIZING,can_inputs,can_outputs,corefunction) = is_k_canalizing_return_inputs_outputs_corefunction(f,n,n)
    return (n,len(can_inputs),can_inputs,can_outputs,corefunction)
   
def get_canalizing_depth_inputs_outputs_corefunction_order(f,variables = []):
    '''obsolete, included for backward compatability, use find_layers(f)'''
    n = int(np.log2(len(f)))
    (NESTED_CANALIZING,can_inputs,can_outputs,corefunction,can_order) = is_k_canalizing_return_inputs_outputs_corefunction_order(f,n,n,variables=variables)
    return (n,len(can_inputs),can_inputs,can_outputs,corefunction,can_order)    


def get_layerstructure_given_canalizing_outputs_and_corefunction(can_outputs,core_polynomial,n=-1):
    '''
    This function computes the canalizing layer structure of a Boolean function 
    given its canalized output values and core function or core polynomial, as defined in [1].
    Two consecutive canalizing variables are in the same (different) layer if 
    they possess the same (different) canalized output value.
    
    Input:
        can_outputs: list of all k canalized output values of a Boolean function f in n variables

        core_polynomial: Boolean core polynomial in all n-k non-conditionally canalizing variables 
        
        n (optional): number of variables of f
        
    Outputs:        
        layerstructure: [k_1,..., k_r] a list of integers describing the number of variables in each canalizing layer of f,
        k_i >= 1, and k_r >= 2 if sum(k_i)==n (i.e, if f is an NCF) unless r = n = 1. 
        
    references:
        [1] Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.

        [2] Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions: Algorithms and applications. Automatica, 146, 110630.
    '''     
    depth = len(can_outputs)
    if depth == 0:
        return []
    if n==-1:
        n = int(np.log2(len(core_polynomial))) + depth
    assert depth!=n-1,"len(can_outputs) == n-1, this is impossible because the last variable is also canalizing in this case."
    if depth == n and n>1: #The last layer of Boolean NCFs has size >=2
        can_outputs[-1] = can_outputs[-2]
    elif is_constant(core_polynomial) and depth>1: #Exceptional case, again last layer here needs to be size >=2
        can_outputs[-1] = can_outputs[-2]
    layerstructure = []
    size_of_layer = 1
    for i in range(1,depth):
        if can_outputs[i]==can_outputs[i-1]:
            size_of_layer+=1
        else:
            layerstructure.append(size_of_layer)
            size_of_layer=1
    layerstructure.append(size_of_layer)
    return layerstructure
    
def get_layerstructure_of_an_NCF_given_its_Hamming_weight(n,w):
    '''
    This function computes the canalizing layer structure of an NCF 
    with a given number of variables and a given Hamming weight. 
    There is a bijection between the Hamming weight 
    (assuming w is equivalent to 2^n-w) and the canalizing layer structure of an NCF.
    
    Input:
        n: number of variables of the NCF,
        
        w: odd Hamming weight of the NCF, i.e., the number of 1s in the NCF in 2^n-vector form
    Outputs:
        r: number of layers of the NCF
        
        layerstructure_NCF: [k_1,..., k_r] a list of integers describing the number of variables in each layer of the NCF,
        k_i >= 1, and k_r >= 2 unless r = n = 1. 
        
    references:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    '''         
    if w==1:
        r=1
        layerstructure_NCF=[n]
    else:
        assert type(w) == int or type(w) == np.int64, 'Hamming weight must be an integer'
        assert 1<=w<=2**n-1, 'Hamming weight w must satisfy, 1 <= w <= 2^n-1'
        assert w%2==1, 'Hamming weight must be an odd integer since all NCFs have an odd Hamming weight.'
        w_bin=dec2bin(w,n)
        current_el=w_bin[0]
        layerstructure_NCF=[1]
        for el in w_bin[1:-1]:
            if el==current_el:
                layerstructure_NCF[-1]+=1
            else:
                layerstructure_NCF.append(1)
                current_el=el
        layerstructure_NCF[-1]+=1
        r=len(layerstructure_NCF)
    return (r,layerstructure_NCF)


## 3) Methods to randomly generate Boolean functions (uniform distribution) and Boolean networks
#Shunting-yard algorithm to evaluate expression
operators = {
    "or": 1,
    "and": 2,
    "not": 3,
    "(": 18,
    ")": 18
}

def eval_expr(expr, x):
    op_stack = []
    val_stack = []
    prevToken = ""
    for token in expr.split(' '):
        if token == '':
            continue

        if token.isdigit():
            val_stack.append(int(token))
        elif not token in operators:
            val = x[int(token[2:-1])]
            val_stack.append(val)
        elif token == '(':
            op_stack.append(token)
        elif token == ')':
            while op_stack[-1] != '(':
                apply_first_op(op_stack, val_stack)

            op_stack.pop()
        else:
            while len(op_stack) > 0 and op_stack[-1] != "(" and get_precedence(op_stack[-1]) >= get_precedence(token):
                apply_first_op(op_stack, val_stack)
            op_stack.append(token)
        
        prevToken = token

    while len(op_stack) > 0:
        apply_first_op(op_stack, val_stack)

    return val_stack[0]

#Helper functions to eval_expr()
def apply_first_op(op_stack, val_stack):
    assert len(op_stack) > 0
    operator = op_stack.pop()
    if operator == "not":
        val_stack.append(int(not val_stack.pop()))
        return

    val1 = val_stack.pop()
    val2 = val_stack.pop()
    outval = 0
    outval = apply_operator(operator, val1, val2)
    val_stack.append(outval)

def get_precedence(operator):
    return operators[operator]
def apply_operator(operator, val1, val2):
    if operator == "not":
        return int(not val1)
    elif operator == "and":
        return int(val1 and val2)
    elif operator == "or":
        return int(val1 or val2)
    else:
        print("err, unrecognized operator: ", operator)

def f_from_expression(expr):
    '''
    This function extracts a Boolean function from a string. 
    
    Input: a text string containing an evaluable expression
        
    Outputs:
        f: the right-hand side of a Boolean function (an array of length 2**n where n is the number of inputs)
    
        var: a list of length n, describing the names and order of variables of f.
        The order is determined by the first occurence of each variable in the input string.
    
    Examples:
        #an and-not function
        print(f_from_expression('A AND NOT B'))
        >> ([0, 0, 1, 0], ['A', 'B'])
        
        #a threshold function
        print(f_from_expression('x1 + x2 + x3 > 1'))
        >> ([0, 0, 0, 1, 0, 1, 1, 1], ['x1', 'x2', 'x3'])
        
        #a parity function
        print(f_from_expression('(x1 + x2 + x3) % 2 == 0'))
        >> ([1, 0, 0, 1, 0, 1, 1, 0], ['x1', 'x2', 'x3'])
    '''    
    expr = expr.replace('(',' ( ').replace(')',' ) ')
    expr_split = expr.split(' ')
    var = []
    dict_var = dict()
    n_var = 0
    for i,el in enumerate(expr_split):
        if el not in ['',' ','(',')','and','or','not','AND','OR','NOT','&','|','~','+','-','*','%','>','>=','==','<=','<'] and not el.isdigit():
            try:
                new_var = dict_var[el]
            except KeyError:
                new_var = 'x[%i]' % n_var
                dict_var.update({el:new_var})
                var.append(el)
                n_var += 1
            expr_split[i] = new_var
        elif el in ['AND','OR','NOT']:
            expr_split[i] = el.lower()
    expr = ' '.join(expr_split)
    f = []
    for x in itertools.product([0, 1], repeat = n_var):
        x = list(map(bool,x))
        f.append(int(eval(expr))) #x is used here "implicitly"
    return f,var

def random_function(n,probability_one=0.5):
    '''
    This function generates a random Boolean function in n variables,
    which are not guaranteed to be essential. 
    
    Inputs:
        n: the number of variables
        
        probability_one (default = 0.5): the bias of the Boolean function,
        i.e., the probability of having a 1 (vs a 0) in the Boolean function.
                
    Output: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    '''   
    return np.array(np.random.random(2**n)<probability_one,dtype=int)
    

def random_linear_function(n):
    '''
    This function generates a random linear Boolean function in n variables.
    
    Input: the number of variables
        
    Output: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    '''   
    return f_from_expression('(%s) %% 2 == 1' % (' + '.join(['x%i' % i if random.random()>0.5 else '(1 + x%i)' % i for i in range(n)])))[0]

def random_non_degenerated_function(n,probability_one=0.5):
    '''
    This function generates a random non-degenerated Boolean function in n variables.
    That is, it generates a Boolean function which possesses only essential variables.
    
    Inputs:
        n: the number of variables
        
        probability_one (default = 0.5): the bias of the Boolean function,
        i.e., the probability of having a 1 (vs a 0) in the Boolean function.
        
    Output: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
        
    references:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    '''  
    while True: #works because most functions are non-degenerated
        f = np.array(np.random.random(2**n)<probability_one,dtype=int)#np.random.randint(2, size = 2**n) 
        if not is_degenerated(f):
            return f

def random_degenerated_function(n,probability_one=0.5):
    '''
    This function generates a random degenerated Boolean function in n variables.
    That is, it generates a Boolean function which possesses at least one non-essential variable.
    
    Inputs:
        n: the number of variables
        
        probability_one (default = 0.5): the bias of the Boolean function,
        i.e., the probability of having a 1 (vs a 0) in the Boolean function.
        
    Output: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)

    references:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    '''  
    while True: #works not well because most functions are non-degenerated
        f = np.array(np.random.random(2**n)<probability_one,dtype=int)#np.random.randint(2, size = 2**n) 
        if is_degenerated(f):
            return f

def random_non_canalizing_function(n,probability_one=0.5):
    '''
    This function generates a random non-canalizing Boolean function in n>1 variables.
    A Boolean function f(x_1,...,x_n) is canalizing if it is canalizing in at least one variable.
    A Boolean function f(x_1,...,x_n) is canalizing in x_i if f(x_1,...,x_i=a,...,x_n) = b for some a,b in [0,1] and for all x_1,...,x_{i-1},x_{i+1},...,x_n in [0,1].

    
    Inputs:
        n: the number of variables
        
        probability_one (default = 0.5): the bias of the Boolean function,
        i.e., the probability of having a 1 (vs a 0) in the Boolean function.
        
    Output: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
    references:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    '''  
    assert n>1
    while True: #works because most functions are non-canalizing
        f = np.array(np.random.random(2**n)<probability_one,dtype=int)#np.random.randint(2, size = 2**n) 
        if not is_canalizing(f,n):
            return f

def random_non_canalizing_non_degenerated_function(n,probability_one=0.5):
    '''
    This function generates a random non-canalizing non-degenerated Boolean function in n>1 variables.
    A Boolean function f(x_1,...,x_n) is canalizing if it is canalizing in at least one variable.
    A Boolean function f(x_1,...,x_n) is canalizing in x_i if f(x_1,...,x_i=a,...,x_n) = b for some a,b in [0,1] and for all x_1,...,x_{i-1},x_{i+1},...,x_n in [0,1].
    A Boolean function is degenerated if at least one of its variables never impacts the output (i.e, is non-essential).
    
    Inputs:
        n: the number of variables (n>1 required)
        
        probability_one (default = 0.5): the bias of the Boolean function,
        i.e., the probability of having a 1 (vs a 0) in the Boolean function.
        
    Output: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
    references:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    ''' 
    assert n>1
    while True: #works because most functions are non-canalizing and non-degenerated
        f = np.array(np.random.random(2**n)<probability_one,dtype=int)#np.random.randint(2, size = 2**n) 
        if not is_canalizing(f,n) and not is_degenerated(f):
            return f

def random_k_canalizing(n, k, EXACT_DEPTH_K=False, left_side_of_truth_table=[]):
    '''
    This function generates a random k-canalizing Boolean function in n variables.
    A Boolean function f is k-canalizing if it has at least k conditionally canalizing variables.

    Inputs:
        n: the number of variables of f
        
        k: the number of conditionally canalizing variables of f
        
        EXACT_DEPTH_K (optional):
            If TRUE, the generated function f has exactly k conditionally canalizing variables, i.e. canalizing depth = k.
            If FALSE (default), the generated function f has canalizing depth >= k.
        
        left_side_of_truth_table (optional for speed-up): the left-hand side of the Boolean truth table of size 2^n x n
        
    Output: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
    references:
        [1] He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth. Physica D: Nonlinear Phenomena, 314, 1-8.

        [2] Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions: Algorithms and applications. Automatica, 146, 110630.
    ''' 
    try:
        assert (n-k!=1 or EXACT_DEPTH_K==False)
    except AssertionError:
        print('There are no functions of exact canalizing depth n-1.\nEither set EXACT_DEPTH_K=False or ensure k!=n-1')
        return
    try:
        assert 0<=k and k<=n
    except AssertionError:
        print('Error:\nEnsure 0 <= k <= n.')
        return
    if left_side_of_truth_table==[]: #to reduce run time, this should be calculated once and then passed as argument
        left_side_of_truth_table = list(itertools.product([0, 1], repeat = n))
    num_values = 2**n
    aas = np.random.randint(2, size = k)  # inputs
    bbs = np.random.randint(2, size = k)  # outputs
    can_vars = np.random.choice(n, k, replace = False)
    f = np.zeros(num_values, dtype = int)
    
    if k<n:
        if EXACT_DEPTH_K:
            core_polynomial = random_non_canalizing_non_degenerated_function(n-k)
        else:
            core_polynomial = random_non_degenerated_function(n-k)    
    else:
        core_polynomial = [1-bbs[-1]]
        
    counter_non_canalized_positions = 0
    
    for i in range(num_values):
        for j in range(k):
            if left_side_of_truth_table[i][can_vars[j]] == aas[j]:
                f[i] = bbs[j]
                break
        else:
            f[i] = core_polynomial[counter_non_canalized_positions]
            counter_non_canalized_positions += 1
    return f

def random_k_canalizing_with_specific_layerstructure(n, layerstructure, EXACT_DEPTH_K=False, left_side_of_truth_table=[]):
    '''
    This function generates a random Boolean function in n variables 
    with a specified layer structure for the k (0<k<=n) canalizing variables.
    
    Inputs:
        n: the number of variables of f
        
        layerstructure: [k_1,..., k_r] a list of integers describing the number of variables in each canalizing layer of f,
        k_i >= 1, and k_r >= 2 if sum(k_i)==n  (i.e, if f is an NCF) unless r = n = 1. 
        
        EXACT_DEPTH_K (optional):
            If TRUE, the generated function f has exactly k conditionally canalizing variables, i.e. canalizing depth = k.
            If FALSE (default), the generated function f has canalizing depth >= k.
        
        left_side_of_truth_table (optional for speed-up): the left-hand side of the Boolean truth table of size 2^n x n
        
    Output: f: a Boolean function as a vector, i.e., the right-hand side of a truth table (a list of length 2^n where n is the number of inputs)
    
    references:
        [1] He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth. Physica D: Nonlinear Phenomena, 314, 1-8.

        [2] Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    ''' 
    k=sum(layerstructure) #canalizing depth
    if k==0:
        layerstructure = [0]
    try:
        assert (n-k!=1 or EXACT_DEPTH_K==False)
    except AssertionError:
        print('Error:\nThere are no functions of exact canalizing depth n-1.\nEither set EXACT_DEPTH_K=False or ensure k=sum(layerstructure)!=n.')
        return
    try:
        assert 0<=k and k<=n
    except AssertionError:
        print('Error:\nEnsure 0 <= k=sum(kis) <= n.')
        return
    try:
        assert k<n or layerstructure[-1]>1 or n==1
    except AssertionError:
        print('Error:\nThe last layer of an n-canalizing function (NCF) has to have size >=2 for n>1.\nIf k=sum(layerstructure)=n, ensure that layerstructure[-1]>=2.')
        return
    try:
        assert min(layerstructure)>=1
    except AssertionError:
        print('Error:\nThere needs to be at least one variable in each layer, i.e., each element of kis must be >=1.')
        return
    if left_side_of_truth_table==[]: #to decrease run time, this should be calculated once and then passed as argument
        left_side_of_truth_table = list(itertools.product([0, 1], repeat = n))
    num_values = 2**n
    aas = np.random.randint(2, size = k)  # inputs
    b0 = np.random.randint(2)
    bbs = [b0]*layerstructure[0]  # outputs
    for i in range(1,len(layerstructure)):
        if i%2==0:
            bbs.extend([b0]*layerstructure[i])
        else:
            bbs.extend([1-b0]*layerstructure[i])
    can_vars = np.random.choice(n, k, replace = False)
    f = np.zeros(num_values, dtype = int)
    
    if k<n:
        if EXACT_DEPTH_K:
            core_polynomial = random_non_canalizing_non_degenerated_function(n-k)
        else:
            core_polynomial = random_non_degenerated_function(n-k)    
    else:
        core_polynomial = [1-bbs[-1]]
        
    counter_non_canalized_positions = 0
    
    for i in range(num_values):
        for j in range(k):
            if left_side_of_truth_table[i][can_vars[j]] == aas[j]:
                f[i] = bbs[j]
                break
        else:
            f[i] = core_polynomial[counter_non_canalized_positions]
            counter_non_canalized_positions += 1
    return f

def random_adj_matrix(N,ns,NO_SELF_REGULATION=True,STRONGLY_CONNECTED=False): #recursive function definition
    matrix = np.zeros((N, N), dtype = int)
    indices = []
    for i in range(N):
        if NO_SELF_REGULATION:
            indexes = np.random.choice( np.append(np.arange(i),np.arange(i+1,N)), ns[i], replace=False)
        else:
            indexes = np.random.choice(np.arange(N), ns[i],replace=False)         
        indexes = np.sort(indexes)
        indices.append(indexes)
        for index in indexes:
            matrix[i][index] = 1
    if STRONGLY_CONNECTED:
        G = nx.from_numpy_array(matrix, parallel_edges = False, create_using = nx.MultiDiGraph())
        if not nx.is_strongly_connected(G):
            return random_adj_matrix(N,ns,NO_SELF_REGULATION,STRONGLY_CONNECTED)
    return (matrix, indices)

def random_edge_list_old(N,ns,NO_SELF_REGULATION):
    edge_list = []
    for i in range(N):
        #indices = np.random.choice(list(range(0, i)) + list(range(i+1, N)),ns[i],replace=False)
        if NO_SELF_REGULATION:
            indexes = np.random.choice( np.append(np.arange(i),np.arange(i+1,N)), ns[i], replace=False)
        else:
            indexes = np.random.choice(np.arange(N), ns[i],replace=False)   
        edge_list.extend( list(zip(indexes,i*np.ones(ns[i],dtype=int))) )
    return edge_list

def random_edge_list(N,ns,NO_SELF_REGULATION,AT_LEAST_ONE_REGULATOR_PER_GENE=False):
    if AT_LEAST_ONE_REGULATOR_PER_GENE==False:
        edge_list = []
        for i in range(N):
            if NO_SELF_REGULATION:
                indices = np.random.choice( np.append(np.arange(i),np.arange(i+1,N)), ns[i], replace=False)
            else:
                indices = np.random.choice(np.arange(N), ns[i],replace=False)   
            edge_list.extend( list(zip(indices,i*np.ones(ns[i],dtype=int))) )
    else:
        edge_list = []
        outdegree = np.zeros(N,dtype=int)
        sum_ns = sum(ns) #total number of regulations
        for i in range(N):
            if NO_SELF_REGULATION:
                indices = np.random.choice( np.append(np.arange(i),np.arange(i+1,N)), ns[i], replace=False)
            else:
                indices = np.random.choice(np.arange(N), ns[i],replace=False)  
            outdegree[indices] += 1
            edge_list.extend( list(zip(indices,i*np.ones(ns[i],dtype=int))) )
        while min(outdegree)==0:
            index_sink = np.where(outdegree==0)[0][0]
            index_edge = int(random.random()*sum_ns)
            if NO_SELF_REGULATION:
                while edge_list[index_edge][1] == index_sink:
                    index_edge = int(random.random()*sum_ns)
            outdegree[index_sink] += 1 
            outdegree[edge_list[index_edge][0]] -= 1 
            edge_list[index_edge] = (index_sink,edge_list[index_edge][1])
    return edge_list

def get_essential_network(F,I):
    '''
    This function determines the number of essential variables of a Boolean function.
    
    inputs:    
        F: a list of N lists of length 2^(n_i) 
        where N is the number of genes and n_i is the number of regulators per gene.
        
        I: a list of N lists of length n_i where N is the number of genes, n_i is the number of regulators per gene. 
        I describes the regulators of each gene (as indices 0, 1, ..., n-1), 
        i.e., the wiring diagram of the Boolean network F (the adjacency matrix of a directed graph)
            
    output:
        F_essential,: a list of N lists of length 2^(m_i), 
        where N is the number of genes and m_i<=n_i is the number of ESSENTIAL regulators per gene.
        
        I_essential: a list of N lists of length m_i 
        where N is the number of genes, m_i<=n_i is the number of ESSENTIAL regulators per gene. 
        I_essential describes the ESSENTIAL regulators of each gene (as indices 0, 1, ..., n-1), 
        i.e., the ESSENTIAL wiring diagram of the Boolean network F (the adjacency matrix of a directed graph)
    '''        
    import itertools
    F_essential = []
    I_essential = []
    for f,regulators in zip(F,I):
        if len(f)==0: #happens if the actual degree of f was too large for it to be loaded
            F_essential.append(f)
            I_essential.append(regulators)
            continue
        elif sum(f)==0:
            F_essential.append(np.array([0]))
            I_essential.append(np.array([],dtype=int))
            continue
        elif sum(f)==len(f):
            F_essential.append(np.array([1]))
            I_essential.append(np.array([],dtype=int))
            continue
        essential_variables = np.array(get_essential_variables(f))
        n = len(regulators)
        non_essential_variables = np.array(list(set(list(range(n))) - set(essential_variables)))
        if len(non_essential_variables)==0:
            F_essential.append(f)
            I_essential.append(regulators)
        else:
            left_side_of_truth_table = np.array(list(itertools.product([0, 1], repeat=n)))
            F_essential.append((np.array(f)[np.sum(left_side_of_truth_table[:,non_essential_variables],1) == 0]))
            I_essential.append((np.array(regulators)[essential_variables]))
    return F_essential,I_essential

def random_BN(N, n = 2, k = 0, STRONGLY_CONNECTED = True, indegree_distribution = 'constant', list_x=[], kis = None, EXACT_DEPTH=False,NO_SELF_REGULATION=True,LINEAR=False,edges_wiring_diagram = None, bias = 0.5):    
    #need to also accept vectors for k and kis
    if indegree_distribution in ['constant','dirac','delta']:
        if type(n) in [list,np.array]:
            try:
                assert type(n) in [list,np.array]
                assert np.all([type(el) in [int,np.int64] for el in n])
                assert len(n)==N
                assert min(n)>=1
                assert max(n)<=N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nTo use a user-defined in-degree vector, ensure that n is an N-dimensional vector where each element of n is an integer between 1 and N.')
                return
        else:
            try:
                assert type(n) in [int,np.int64]
                assert n>=1
                assert n<=N
                ns = np.ones(N,dtype=int)*n
            except AssertionError:
                print('Error: n must be a single integer (or N-dimensional vector of integers) between 1 and N when using a constant degree distribution.')
                return            
    elif indegree_distribution == 'uniform':
        if type(n) in [list,np.array]:
            try:
                assert type(n) in [list,np.array]
                assert np.all([type(el) in [int,np.int64] for el in n])
                assert len(n)==N
                assert min(n)>=1
                assert max(n)<=N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nEnsure that n is an N-dimensional vector where each element of n is an integer between 1 and N representing the upper bound of a uniform degree distribution (lower bound==1).')
                return
        else:
            try:
                assert type(n) in [int,np.int64]
                assert n>=1
                assert n<=N
                ns = np.ones(N,dtype=int)*n
            except AssertionError:
                print('Error: n must be a single integer (or N-dimensional vector of integers) between 1 and N representing the upper bound of a uniform degree distribution (lower bound==1).')
                return
    elif indegree_distribution == 'poisson':
        if type(n) in [list,np.array]:
            try:
                assert type(n) in [list,np.array]
                assert np.all([type(el) in [int,np.int64,float,np.float64] for el in n])
                assert len(n)==N
                assert min(n)>=1
                assert max(n)<=N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nEnsure that n is an N-dimensional vector where each element of n is > 0 and represents the Poisson parameter.')
                return
        else:
            try:
                assert type(n) in [int,np.int64,float,np.float64]
                assert n>=1
                assert n<=N
                ns = np.ones(N,dtype=int)*n
            except AssertionError:
                print('Error: n must be a single number (or N-dimensional vector) > 0 representing the Poisson parameter.')
                return            
    else:
        print('None of the predefined indegree distributions were chosen.\nTo use a user-defined in-degree vector, use the input n to submit an N-dimensional vector where each element of n must be between 1 and N')
        return

    if kis==None:
        if type(k) in [int,np.int64]:
            try:
                assert k>=0
                assert k<=N
                max_k = k
            except AssertionError:
                print('Error: k must be an integer between 0 and N.')
                return
        elif type(k) in [list,np.array]:
            try:
                assert len(k)==N
                assert np.all([type(el) in [int,np.int64] for el in k])
                max_k = max(k)
                assert min(k)>=0
                assert max_k<=N
            except AssertionError:
                print('Error: A vector k was submitted.\nTo use a user-defined vector k, ensure that k is an N-dimensional vector where each element of k is an integer between 0 and N.')
                return
        else:
            print('Error: Wrong input format for k.\nk must be a single integer (or N-dimensional vector of integers) between 0 and N')
            return
    else: #kis provided
        if np.all([type(el) in [int,np.int64] for el in kis]):
            try:
                assert np.all([type(el) in [int,np.int64] for el in kis])
                assert np.all([el>=1 for el in kis])
                max_k = sum(kis)
                assert max_k<=N #could be more restrictive but if n is also a vector this is tedious to program
            except AssertionError:
                print('Error: the layer structure kis must be a vector of positive integers with 0<= k=sum(kis) <= n.')
                return
        elif np.all([type(el) in [list,np.array] for el in kis]):
            try:
                assert len(kis)==N
                assert type(kis[0][0]) in [int,np.int64]
                max_k = max([sum(el) for el in kis])
                assert min([min(el) for el in kis])>=0
                assert max_k<=N
            except AssertionError:
                print('Error: A vector of kis was submitted.\nTo use a user-defined vector of kis, ensure that kis is an N-dimensional vector where each element represents a layer structure and is a vector of positive integers with 1<= k=sum(kis[i]) <= n.')
                return
        else:
            print('Error: Wrong input format for kis.\nkis must be a single vector (or N-dimensional vector of layer structures) where each the sum of each element must be between 0 and N.')
            return

    if edges_wiring_diagram is None:
        while True: # Keep generating until we have a strongly connected graph
            # create a new degree vector
            # Can be constant, random or a statistical distribution
            if indegree_distribution == 'uniform':            
                ns = 1 + np.random.randint(n-1, size = N)
            elif indegree_distribution == 'poisson':
                ns = np.random.poisson(lam = n , size = N)
                ns[ns==0] = 1
                ns[ns>N-int(NO_SELF_REGULATION)] = N-int(NO_SELF_REGULATION)
                
            #A, I = random_adj_matrix(N, ns, NO_SELF_REGULATION)
            edges_wiring_diagram = random_edge_list(N,ns,NO_SELF_REGULATION)
            
            # If we care whether the graph is strongly connected (STRONGLY_CONNECTED),
            # we will check the connectedness of the graph using nx.is_strongly_connected(G)
            if STRONGLY_CONNECTED:
                #G = nx.from_numpy_array(A, parallel_edges = False, create_using = nx.MultiDiGraph())
                G = nx.from_edgelist(edges_wiring_diagram, create_using = nx.MultiDiGraph())
                strongly_connected = nx.is_strongly_connected(G)
                if not strongly_connected:
                    continue
            break
    else:
        try:
            assert len(set(np.array(edges_wiring_diagram).flatten())) == N
        except AssertionError:
            print("number of nodes provided in edges_wiring_diagram != N")
            return
        ns = np.zeros(N,dtype=int)
        for target in np.array(edges_wiring_diagram)[:,1]:
            ns[target] += 1
        
    max_n = max(ns)
    if max_k>0 and list_x==[] or len(list_x)<max_n: #list_x probably not correct, or fully generated
        #create list_x
        list_x = [[[0],[1]]]
        list_x.extend([list(itertools.product([0, 1], repeat = nn)) for nn in range(2, max_n+1)])

        
    F = []
    for i in range(N):
        if LINEAR:
            F.append(random_linear_function(ns[i]))
        if k>0 and kis==None:
            if type(k) in [int,np.int64]:
                F.append(random_k_canalizing(ns[i], min(k,ns[i]), EXACT_DEPTH_K = EXACT_DEPTH, x=list_x[ns[i]-1]))
            else:
                F.append(random_k_canalizing(ns[i], min(k[i],ns[i]), EXACT_DEPTH_K = EXACT_DEPTH, x=list_x[ns[i]-1]))
        elif kis!=None: #value of k is ignored if a layer structure is provided
            if np.all([type(el) in [int,np.int64] for el in kis]):
                F.append(random_k_canalizing_with_specific_layerstructure(ns[i], kis, EXACT_DEPTH_K = EXACT_DEPTH, x=list_x[ns[i]-1]))                
            else:
                F.append(random_k_canalizing_with_specific_layerstructure(ns[i], kis[i], EXACT_DEPTH_K = EXACT_DEPTH, x=list_x[ns[i]-1]))                                
        else:
            if EXACT_DEPTH==True: #i.e. if k==0
                F.append(random_non_canalizing_non_degenerated_function(ns[i],bias))                   
            else:
                F.append(random_non_degenerated_function(ns[i],bias))
    
    I = [[] for _ in range(N)]
    for edge in edges_wiring_diagram:
        I[edge[1]].append(edge[0])
    for i in range(N):
        I[i] = np.sort(I[i])
    
    return F, I, ns

def get_perturbed_network(F,I,ns,control_target,control_source,type_of_control = 0,left_side_of_truth_table = []):
    F_new = [f for f in F]
    I_new = [i for i in I]
    
    if left_side_of_truth_table==[]:
        left_side_of_truth_table = np.array(list(itertools.product([0,1],repeat=ns[control_target])))
    
    try:
        index = list(I[control_target]).index(control_source)
        F_new[control_target] = F_new[control_target][left_side_of_truth_table[:,index]==type_of_control]
        dummy = list(I_new[control_target])
        dummy.remove(control_source)
        I_new[control_target] = np.array(dummy)
    except ValueError:
        print('source not regulating target')
    
    ns = list(map(len,I_new))
    return F_new,I_new,ns


# E = random_edge_list(N,ns,NO_SELF_REGULATION)

# # If we care whether the graph is strongly connected (STRONGLY_CONNECTED),
# # we will check the connectedness of the graph using nx.is_strongly_connected(G)
# if STRONGLY_CONNECTED:
#     #G = nx.from_numpy_array(A, parallel_edges = False, create_using = nx.MultiDiGraph())
#     G = nx.from_edgelist(E, create_using = nx.MultiDiGraph())
#     strongly_connected = nx.is_strongly_connected(G)



## 4) Enumeration methods
def nr_non_canalizing_by_weight_exact(n):
    assert n<=4
    
    nr_values = 2**n
    nr_boolean_fcts = 2**(nr_values)

    F = np.array(list(itertools.product([0, 1], repeat=nr_values))).T
    
    ws_can,ws_non_can = [],[]
    for i in range(nr_boolean_fcts):
        if is_canalizing(F[:,i],n):
            ws_can.append(sum(F[:,i]))
        else:
            ws_non_can.append(sum(F[:,i]))
    return [ws_non_can,ws_can]

def nr_non_canalizing_by_weight_simulation(n,nsim=10000):
    nr_values = 2**n
    
    ws_non_can,ws_can = [],[]
    
    for i in range(nsim):
        F = np.random.randint(2, size = nr_values)
        if is_canalizing(F,n):
            ws_can.append(sum(F))
        else:
            ws_non_can.append(sum(F))
    return [ws_non_can,ws_can]

def stratify_Boolean_fcts_by_canalization_ns(n,nsim=10000):
    nr_values = 2**n  
    depths,layerstructures = [],[]
    for i in range(nsim):
        f = np.random.random(nr_values)>0.5
        (depth_f,number_layers_f,can_inputs_f,can_outputs_f,corepolynomial_f,can_order_f) = find_layers(f)
        depths.append(depth_f)
        layerstructure = get_layerstructure_given_canalizing_outputs_and_corefunction(can_outputs_f,corepolynomial_f,n)
        layerstructures.append(layerstructure)
    return (depths,layerstructures)


## 5) Analysis methods
def get_constant_nodes(I,degree,N):
    return np.array([i for i in range(N) if degree[i]==1 and I[i][0]==i])

def rnd_edge_deletion(F,I,N,degree,nsim=100,left_sides_of_truth_tables=[]):
    res = []
    for i in range(nsim):
        rnd_var = np.random.randint(N)
        rnd_input = np.random.randint(degree[rnd_var])
        res.append( sum(F[rnd_var][np.array(left_sides_of_truth_tables[degree[rnd_var]-1])[:,rnd_input]==0]!=F[rnd_var][np.array(left_sides_of_truth_tables[degree[rnd_var]-1])[:,rnd_input]==1]) )
    return res

def update(F, I, N, X):
    Fx = np.zeros(N, dtype = int)
    for i in range(N):
        Fx[i] = F[i][bin2dec(X[I[i]])]
    return Fx    

def average_sensitivity_old_wrong(F,nsim=10000): 
    #equals Derrida value D(F,1) if all n update rules in F are chosen from the same sampling space
    if type(F)==list:
        F = np.array(F)
    n = int(np.log2(len(F)))
    num_values = 2**n
    X = np.random.randint(num_values,size=nsim)
    add = 2**np.random.randint(n, size=nsim)
    Y = (X + add )%num_values
    return sum(np.bitwise_xor(F[X],F[Y]))*1./nsim

def average_sensitivity(F,nsim=10000,EXACT=False,NORMALIZED=True): 
    #equals Derrida value D(F,1) if all n update rules in F are chosen from the same sampling space
    if type(F)==list:
        F = np.array(F)
    n = int(np.log2(len(F)))
    num_values = 2**n
    s = 0
    if EXACT:
        left_side_of_truth_table = list(map(np.array,list(itertools.product([0, 1], repeat = n))))
        for ii,X in enumerate(left_side_of_truth_table):
            for i in range(n):
                Y=X.copy()
                Y[i] = 1-X[i]
                Ydec = bin2dec(Y)
                s += int(F[ii]!=F[Ydec])
        if NORMALIZED:
            return s/num_values/n
        else:
            return s/num_values
    else:
        for i in range(nsim):
            Xbin = np.random.randint(num_values)
            Y = dec2bin(Xbin,n)
            index = np.random.randint(n)
            Y[index] = 1-Y[index]
            Ybin = bin2dec(Y)
            s += int(F[Xbin]!=F[Ybin])
        if NORMALIZED:
            return s/nsim
        else:
            return n*s/nsim


def absolute_bias(F,n=None):
    if n==None:
        n = int(np.log2(len(F)))
    return abs(sum(F)*1./2**(n-1)-1)

def num_of_attractors(F, I, N, nsim = 500, EXACT = False, left_side_of_truth_table = [], initial_sample_points = []):
    dictF = dict()
    attractors = []
    basin_sizes = []
    attr_dict = dict()
    
    if EXACT and left_side_of_truth_table == []:
        left_side_of_truth_table = list(map(np.array,list(itertools.product([0, 1], repeat = N))))
    
    sampled_points = []
    
    for i in range(nsim if not EXACT else 2**N):
        if initial_sample_points==[]:
            x = np.random.randint(2, size = N) if not EXACT else left_side_of_truth_table[i]
            xbin = bin2dec(x)
            sampled_points.append(xbin)
        else:
            xbin = initial_sample_points[i]
            x = dec2bin(xbin,N)
        queue = [xbin]        
        while True: #as long as we haven't reached an attractor state, keep updating
            try:
                fxbin = dictF[xbin]
            except KeyError:
                fx = update(F, I, N, x)
                fxbin = bin2dec(fx)
                dictF.update({xbin:fxbin})
                x = fx #works only because if we don't know fx now we also won't know F[fx] 
            try: # check if this state is a known attractor
                index_attr = attr_dict[fxbin] #returns a KeyError if we haven't identified fxbin as an attractor before
                basin_sizes[index_attr] += 1
                break
            except KeyError:
                try: #check if fxbin is part of a new attractor
                    index = queue.index(fxbin) #returns a ValueError if fxbin is not part of queue
                    #new attractor
                    attr_dict.update( list(zip( queue[index:] , [len(attractors)]*(len(queue)-index) )) )
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(fxbin)
            xbin = fxbin
    return (attractors, len(attractors), basin_sizes, attr_dict, initial_sample_points if initial_sample_points!=[] else sampled_points)

def num_of_attractors_v2(F, I, N, nb = 500, initial_sample_points=[]): #should be faster if we have longer average path lengths, not significantly though
    dictF = dict()
    attractors = []
    basin_sizes = []
    attr_dict = dict()
    
    sampled_points = []
    
    for i in range(nb):
        if initial_sample_points==[]:
            x = np.random.randint(2, size = N)
            xbin = bin2dec(x)
            sampled_points.append(xbin)
        else:
            xbin = initial_sample_points[i]
            x = dec2bin(xbin,N)
        xbin = bin2dec(x)
        queue = [xbin]        
        while True: #as long as we haven't reached an attractor state, keep updating
            try:
                fxbin = dictF[xbin]
            except KeyError:
                fx = update(F, I, N, x)
                fxbin = bin2dec(fx)
                dictF.update({xbin:fxbin})
                x = fx #works only because if we don't know fx now we also won't know F[fx] 
            try: # check if this state has a known attractor
                index_attr = attr_dict[fxbin] #returns a KeyError if we haven't identified fxbin as an attractor before
                basin_sizes[index_attr] += 1
                attr_dict.update( list(zip( queue , [index_attr]*len(queue) )) )
                break
            except KeyError:
                try: #check if fxbin is part of a new attractor
                    index = queue.index(fxbin) #returns a ValueError if fxbin is not part of queue
                    #new attractor
                    attr_dict.update( list(zip( queue[index:] , [len(attractors)]*(len(queue)-index) )) )
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(fxbin)
            xbin = fxbin
    return (attractors, len(attractors), basin_sizes, attr_dict, initial_sample_points if initial_sample_points!=[] else sampled_points)

def num_of_attractors_exact(F, I, N,left_side_of_truth_table = []):
    dictF = dict()
    attractors = []
    basin_sizes = []
    attractor_dict = dict()
    
    if left_side_of_truth_table == []:
        left_side_of_truth_table = list(map(np.array,list(itertools.product([0, 1], repeat = N))))
    
    b_for_bin2decs = [np.array([2**i for i in range(NN)])[::-1] for NN in range(N+1)]
    
    degrees = list(map(len,I))
    
    for xbin,x in enumerate(left_side_of_truth_table):
        queue = [xbin]
        while True: #as long as we haven't reached an attractor state, keep updating
            try:
                fxbin = dictF[xbin]
            except KeyError:
                fx = []
                for i in range(N):
                    #fx.append(F[i][sum([x[I[i]][degrees[i]-j-1]*b_for_bin2dec[j] for j in range(degrees[i])])])
                    #fx.append(F[i][sum([x[I[i]][j]*b_for_bin2dec[j] for j in range(degrees[i])])])
                    #fx.append(F[i][np.dot(x[I[i]] , b_for_bin2dec[N-degrees[i]:])])
                    fx.append(F[i][np.dot(x[I[i]] , b_for_bin2decs[degrees[i]])])
                    
                #fx = update(F,I,N,x)
                fxbin = np.dot(fx,b_for_bin2decs[-1])#sum([fx[N-i-1]*b_for_bin2dec[i] for i in range(N)])
                dictF.update({xbin:fxbin})
                x = np.array(fx) #works only because if we don't know fx now we also won't know F[fx] 
            try: # check if this state has a known attractor
                index_attr = attractor_dict[fxbin] #returns a KeyError if we haven't identified fxbin as an attractor before
                basin_sizes[index_attr] += 1
                attractor_dict.update( list(zip( queue , [index_attr]*len(queue) )) )
                break
            except KeyError:
                try: #check if fxbin is part of a new attractor
                    index = queue.index(fxbin) #returns a ValueError if fxbin is not part of queue
                    #new attractor
                    attractor_dict.update( list(zip( queue , [len(attractors)]*len(queue) )) )
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(fxbin)
            xbin = fxbin
    return (attractors, len(attractors), basin_sizes, attractor_dict)

def num_of_attractors_exact_fast(F, I, N,left_side_of_truth_table = None):    
    if left_side_of_truth_table is None:
        left_side_of_truth_table = np.array(list(map(np.array,list(itertools.product([0, 1], repeat = N)))))
    
    b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
    degrees = list(map(len,I))
    
    state_space = np.zeros((2**N,N),dtype=int)
    for i in range(N):
        for j,x in enumerate(itertools.product([0, 1], repeat = degrees[i])):
            if F[i][j]:
                state_space[np.all(left_side_of_truth_table[:,I[i]]==np.array(x),1),i] = 1
    dictF = dict(zip(list(range(2**N)),np.dot(state_space,b_for_bin2dec)))
    
    attractors = []
    basin_sizes = []
    attractor_dict = dict()    
    for xbin in range(2**N):
        queue = [xbin]
        while True: #as long as we haven't reached an attractor state, keep updating
            fxbin = dictF[xbin]
            try: # check if this state has a known attractor
                index_attr = attractor_dict[fxbin] #returns a KeyError if we haven't identified fxbin as an attractor before
                basin_sizes[index_attr] += 1
                attractor_dict.update( list(zip( queue , [index_attr]*len(queue) )) )
                break
            except KeyError:
                try: #check if fxbin is part of a new attractor
                    index = queue.index(fxbin) #returns a ValueError if fxbin is not part of queue
                    #new attractor
                    attractor_dict.update( list(zip( queue , [len(attractors)]*len(queue) )) )
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(fxbin)
            xbin = fxbin
    return (attractors, len(attractors), basin_sizes, attractor_dict)
        
def basin_size_largest(basin_sizes):
    return max(basin_sizes)*1./sum(basin_sizes)

def entropy(basin_sizes):
    total = sum(basin_sizes)
    return sum([-np.log(el)*el for el in [size*1./total for size in basin_sizes]])

def d(x, y):            
    return sum(x!=y) #works as long as x,y are numpy arrays

def derrida_value(F, I, N, m, nsim = 500):
    total = 0
    for i in range(nsim):
        X = np.random.randint(2, size = N)
        ones = np.zeros(N, dtype = int)
        ones[np.random.choice(N, m, replace = False)] = 1
        Y = np.bitwise_xor(X, ones)
        total += d(update(F, I, N, X), update(F, I, N, Y))
    return total*1./nsim

def get_robustness_from_attractor_dict_exact(attractor_dict,N,n_attractors,left_side_of_truth_table):
    '''computes the proportion of neighbors in the Boolean hypercube
    who, following synchronous update, transition to the same attractor,
    
    Note: attr_dict must contain a fully sampled state space'''
    if n_attractors==1:
        return 1
    b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
    count_of_neighbors_who_transition_to_same_attractor = 0
    for xbin,x in enumerate(left_side_of_truth_table):
        for i in range(N):
            if x[i]==0:
                ybin = xbin + b_for_bin2dec[i]
            else:
                continue
            if attractor_dict[xbin] == attractor_dict[ybin]:
                count_of_neighbors_who_transition_to_same_attractor += 1
    return count_of_neighbors_who_transition_to_same_attractor/2**(N-1)/N

def get_robustness_and_attractors_simulation(F, I, N, number_different_IC = 500): 
    '''computes by sampling the attractor landscape as well as
    the proportion of neighbors in the Boolean hypercube
    who, following synchronous update, transition to the same attractor'''
    dictF = dict()
    attractors = []
    basin_sizes = []
    attractor_dict = dict()

    degrees = list(map(len,I))
    
    b_for_bin2decs = [np.array([2**i for i in range(NN)])[::-1] for NN in range(max(degrees)+1)]
    b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
    
    robustness_approximation = 0
    
    for i in range(number_different_IC):
        index_attractors = []
        for j in range(2):
            if j==0:
                x = np.random.randint(2, size = N)
                xbin = np.dot(x,b_for_bin2dec)
                x_old = x.copy()
            else:
                x = x_old
                random_flipped_bit = np.random.choice(N)
                x[random_flipped_bit] = 1 - x[random_flipped_bit]
                xbin = np.dot(x,b_for_bin2dec)                
            queue = [xbin]
            while True: #as long as we haven't reached an attractor state, keep updating
                try:
                    fxbin = dictF[xbin]
                except KeyError:
                    fx = []
                    for jj in range(N):
                        #fx.append(F[i][sum([x[I[i]][degrees[i]-j-1]*b_for_bin2dec[j] for j in range(degrees[i])])])
                        #fx.append(F[i][sum([x[I[i]][j]*b_for_bin2dec[j] for j in range(degrees[i])])])
                        #fx.append(F[i][np.dot(x[I[i]] , b_for_bin2dec[N-degrees[i]:])])
                        fx.append(F[jj][np.dot(x[I[jj]] , b_for_bin2decs[degrees[jj]])])
                        
                    #fx = update(F,I,N,x)
                    fxbin = np.dot(fx,b_for_bin2dec)#sum([fx[N-i-1]*b_for_bin2dec[i] for i in range(N)])
                    dictF.update({xbin:fxbin})
                    x = np.array(fx) #works only because if we don't know fx now we also won't know F[fx] 
                try: # check if this state has a known attractor
                    index_attr = attractor_dict[fxbin] #returns a KeyError if we haven't identified fxbin as an attractor before
                    basin_sizes[index_attr] += 1
                    attractor_dict.update( list(zip( queue , [index_attr]*len(queue) )) )
                    break
                except KeyError:
                    try: #check if fxbin is part of a new attractor
                        index = queue.index(fxbin) #returns a ValueError if fxbin is not part of queue
                        #new attractor
                        index_attr = len(attractors)
                        attractor_dict.update( list(zip( queue , [index_attr]*len(queue) )) )
                        attractors.append(queue[index:])
                        basin_sizes.append(1)
                        break
                    except ValueError:
                        pass
                queue.append(fxbin)
                xbin = fxbin
            index_attractors.append(index_attr)
        if index_attractors[0] == index_attractors[1]:
            robustness_approximation += 1
        
    return (attractors, len(attractors), basin_sizes, robustness_approximation/number_different_IC)





def adjacency_matrix(I,constants=[],IGNORE_SELFLOOPS=False,IGNORE_CONSTANTS=True):
    n = len(I)
    n_constants = len(constants)
    if IGNORE_CONSTANTS:
        m = np.zeros((n-n_constants,n-n_constants),dtype=int)
        for i,regulators in enumerate(I):
            for j in regulators:
                if j<n-n_constants and (IGNORE_SELFLOOPS==False or i!=j):
                    m[j,i] = 1
        return m
    else:
        return adjacency_matrix(I,[],IGNORE_CONSTANTS=True)

def is_monotonic(F,GET_DETAILS=False):
    n=int(np.log2(len(F)))
    F = np.array(F)
    monotonic = []
    for i in range(n):
        dummy_add=(2**(n-1-i))
        dummy=np.arange(2**n)%(2**(n-i))//dummy_add
        diff = F[dummy==1]-F[dummy==0]
        min_diff = min(diff)
        max_diff = max(diff)
        if min_diff==0 and max_diff==0:
            monotonic.append('not essential')
        elif min_diff==-1 and max_diff==1:
            monotonic.append('not monotonic')
        elif min_diff>=0 and max_diff==1:
            monotonic.append('increasing')            
        elif min_diff==-1 and max_diff<=0:
            monotonic.append('decreasing')   
    if GET_DETAILS:
        return ('not essential' not in monotonic,monotonic)
    else:
        return 'not essential' not in monotonic

def get_ffls(adjacency_matrix,F=None,I=None):
    n = len(adjacency_matrix)
    ffls = []
    for i in range(n):
        for j in range(n):
            if i==j or adjacency_matrix[i,j]==0:
                continue
            #else: find other nodes regulated by i 
            for k in range(n):
                if j==k or i==k or adjacency_matrix[i,k]==0:
                    continue
                #else: check if that other node k also regulates i
                if adjacency_matrix[k,j]!=0: #found FFL
                    ffls.append([i,k,j])
    if F==None or I==None:
        return ffls
    dict_monotonic = dict()
    types = []
    for [j,k,i] in ffls:
        try:
            monotonic = dict_monotonic[k]
        except KeyError:
            monotonic = is_monotonic(F[k],True)[1]
            dict_monotonic.update({k : monotonic})
        indirect1 = monotonic[list(I[k]).index(j)]
        try:
            monotonic = dict_monotonic[i]
        except KeyError:
            monotonic = is_monotonic(F[i],True)[1]
            dict_monotonic.update({i : monotonic})
        indirect2 = monotonic[list(I[i]).index(k)]
        direct = monotonic[list(I[i]).index(j)]
        types.append([direct,indirect1,indirect2])
    return (ffls,types)

def get_ffls_from_I(I,types_I=None):
    all_tfs = list(range(len(I)))#list(set.union(*[] + [set(el) for el in I]))
    n_tfs = len(all_tfs)
    all_tfs_dict = dict(zip(all_tfs,list(range(n_tfs))))
    I_inv = [[] for _ in all_tfs]
    for target,el in enumerate(I):
        for regulator in el:
            I_inv[all_tfs_dict[regulator]].append(target)
    ffls = []
    types = []
    for i in range(n_tfs): #master regulators
        for j in range(n_tfs):
            if i==j or all_tfs[j] not in I_inv[i]:
                continue
            #else: find common nodes regulated by both
            common_targets = list(set(I_inv[i]) & set(I_inv[j]))
            for k in common_targets:
                if all_tfs[j]==k or all_tfs[i]==k:
                    continue
                ffls.append([i,j,k])
                if types_I != None:
                    direct = types_I[k][I[k].index(all_tfs[i])]
                    indirect1 = types_I[all_tfs[j]][I[all_tfs[j]].index(all_tfs[i])]
                    indirect2 = types_I[k][I[k].index(all_tfs[j])]
                    types.append([direct,indirect1,indirect2])
    if types_I != None:
        return (ffls,types)
    else:
        return ffls
    
def get_ffl_type_number(types_vector):
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return -1 if 'not essential' not in types_vector else -2
    else:
        dummy = np.array([1 if el=='increasing' else 0 for el in types_vector])
        nr_type = np.dot(dummy,2**np.arange(len(types_vector)))
    return nr_type

def is_ffl_coherent(types_vector):
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return np.nan
    else:
        dummy = np.array([1 if el=='increasing' else 0 for el in types_vector])
        COHERENT = sum(dummy)%2==1
    return COHERENT

def generate_networkx_graph(I,constants,variables):
    names = list(variables)+list(constants)
    G=nx.DiGraph()
    G.add_nodes_from(names)
    G.add_edges_from([(names[I[i][j]],names[i]) for i in range(len(variables)) for j in range(len(I[i]))])
    return G

def generate_networkx_graph_from_edges(I,n_variables):
    edges = []
    for j,regulators in enumerate(I):
        if j>=n_variables: #exclude constant self-loops
            break
        for i in regulators:
            edges.append((i,j))
    return nx.DiGraph(edges)

def simple_cycles(G,max_len=4):
    '''
    from https://networkx.org/documentation/networkx-1.9/_modules/networkx/algorithms/cycles.html#simple_cycles
    
    This function 
    '''
    def _unblock(thisnode, blocked, B):
        stack = set([thisnode])
        while stack:
            node = stack.pop()
            if node in blocked:
                blocked.remove(node)
                stack.update(B[node])
                B[node].clear()
    
    # Johnson's algorithm requires some ordering of the nodes.
    # We assign the arbitrary ordering given by the strongly connected comps
    # There is no need to track the ordering as each node removed as processed.
    # Also we save the actual graph so we can mutate it. We only take the
    # edges because we do not want to copy edge and node attributes here.
    subG = type(G)(G.edges())
    sccs = [scc for scc in nx.strongly_connected_components(subG)
            if len(scc) > 1]
    
    # Johnson's algorithm exclude self cycle edges like (v, v)
    # To be backward compatible, we record those cycles in advance
    # and then remove from subG
    for v in subG:
        if subG.has_edge(v, v):
            yield [v]
            subG.remove_edge(v, v)
    
    while sccs:
        scc = sccs.pop()
        sccG = subG.subgraph(scc)
        # order of scc determines ordering of nodes
        startnode = scc.pop()
        # Processing node runs "circuit" routine from recursive version
        path = [startnode]
        len_path = 1
        blocked = set()  # vertex: blocked from search?
        closed = set()   # nodes involved in a cycle
        blocked.add(startnode)
        B = defaultdict(set)  # graph portions that yield no elementary circuit
        stack = [(startnode, list(sccG[startnode]))]  # sccG gives comp nbrs
        while stack:
            thisnode, nbrs = stack[-1]
            if nbrs and len_path<=max_len:
                nextnode = nbrs.pop()
                if nextnode == startnode:
                    yield path[:]
                    closed.update(path)
                elif nextnode not in blocked:
                    path.append(nextnode)
                    len_path+=1
                    stack.append((nextnode, list(sccG[nextnode])))
                    closed.discard(nextnode)
                    blocked.add(nextnode)
                    continue
            # done with nextnode... look for more neighbors
            if not nbrs or len_path>max_len:  # no more nbrs
                if thisnode in closed:
                    _unblock(thisnode, blocked, B)
                else:
                    for nbr in sccG[thisnode]:
                        if thisnode not in B[nbr]:
                            B[nbr].add(thisnode)
                stack.pop()
                path.pop()
                len_path-=1
        # done processing this node
        H = subG.subgraph(scc)  # make smaller to avoid work in SCC routine
        sccs.extend(scc for scc in nx.strongly_connected_components(H)
                    if len(scc) > 1)

def get_type_of_loop(loop,F,I):
    n = len(loop)
    dummy = loop[:]
    dummy.append(loop[0])
    res = []
    for i in range(n):
        res.append( is_monotonic(F[dummy[i+1]],True)[1][list(I[dummy[i+1]]).index(dummy[i])] )
    return res

def get_loop_type_number(types_vector):
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return -1 if 'not essential' not in types_vector else -2
    else:
        nr_type = int(np.sum([1 if el=='decreasing' else 0 for el in types_vector]))
    return nr_type

def is_pos_loop(types_vector):
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return np.nan
    else:
        POSITIVE = int(np.sum([1 if el=='decreasing' else 0 for el in types_vector])) % 2 == 0
    return POSITIVE







## 6) Interesting plots
def bias_vs_sens(n,nsim_per_k=1000,avg_sens_sim=1000): 
    assert n>1
    f,ax = plt.subplots()
    colorlist = np.array(list(set(mcolors.BASE_COLORS.keys())-set('w')))
    biass,senss = [],[]
    ks = (list(range(n-1)) + [n])
    for k in ks:
        biass.append([])
        senss.append([])
        for i in range(nsim_per_k):
            F = random_k_canalizing(n,k,True)#random_non_degenerated_function(n)
            biass[-1].append(absolute_bias(F,n))
            senss[-1].append(average_sensitivity(F,avg_sens_sim))
        ax.plot(biass[-1],senss[-1],'o',color=colorlist[k],label=str(k),alpha=0.4)
    ax.legend(loc='best',title='canalizing depth')
    ax.set_xlabel('absolute bias')
    ax.set_ylabel('average sensitivity')
    ax.set_title('Boolean functions in %i essential variables' % n)
    plt.savefig('bias_vs_sens_v01_n%i_nsim%i_avgsenssim%i.pdf' % (n,nsim_per_k,avg_sens_sim))

def bias_vs_sens_v2(n,nsim_per_k=1000,avg_sens_sim=1000): 
    assert n>1
    f,ax = plt.subplots()
    colorlist = np.array(list(set(mcolors.BASE_COLORS.keys())-set('w'))+['orange'])
    biass,senss = [],[]
    ks = (list(range(n-1)) + [n])
    for k in ks:
        biass.append([])
        senss.append([])
        for i in range(nsim_per_k):
            F = random_k_canalizing(n,k,True)#random_non_degenerated_function(n)
            biass[-1].append(absolute_bias(F,n))
            senss[-1].append(average_sensitivity(F,avg_sens_sim))
        which_biases = list(set(biass[-1]))
        which_biases.sort()
        median_sens = [np.median(np.array(senss[-1])[np.array(biass[-1]) == el]) for el in which_biases]
        ax.plot(which_biases,median_sens,'x-',color=colorlist[k],label=str(k),alpha=0.4)
    ax.legend(loc='best',title='canalizing depth')
    ax.set_xlabel('absolute bias')
    ax.set_ylabel('median average sensitivity')
    ax.set_title('Boolean functions in %i essential variables' % n)
    f.savefig('bias_vs_sens_v02_n%i_nsim%i_avgsenssim%i.pdf' % (n,nsim_per_k,avg_sens_sim))
    
    g,ax = plt.subplots()
    ax.boxplot(senss)
    ax.set_xticklabels(ks)
    ax.set_xlabel('canalizing depth')
    ax.set_ylabel('average sensitivity')
    ax.set_title('Boolean functions in %i essential variables' % n)
    g.savefig('average_sensitivity_v02_n%i_nsim%i_avgsenssim%i.pdf' % (n,nsim_per_k,avg_sens_sim))
    
