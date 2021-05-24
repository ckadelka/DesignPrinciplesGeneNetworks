#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import pickle
import os 

folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '27765040'

F, I, degree, var, constants = [],[],[],[],[]
dict_var_and_constants = dict()
count=0
for fname in os.listdir(folder+pmid+'/'):
    if fname.endswith('.txt'):
        var.append(fname.split('.txt')[0])
        dict_var_and_constants.update({var[-1]:count})
        count+=1
        
for fname in os.listdir(folder+pmid+'/'):
    if fname.endswith('.txt'):        
        f = open(folder+pmid+'/'+fname,'r')
        text = f.read().splitlines()
        f.close()
        
        rule = list(map(int,[el.split('\t')[-1] for el in text[1:]]))
        F.append(rule)
        
        regulators = text[0].split('\t')[:-1]

        I.append([])
        for regulator in regulators:
            try:
                index = dict_var_and_constants[regulator]
            except KeyError:
                dict_var_and_constants.update({regulator:count})
                index = count
                count+=1
                constants.append(regulator)
            I[-1].append(index)

g = open(folder+pmid+'_tabular.txt','wb')
pickle.dump([F,I,var,constants],g)
g.close()