#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 

folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '27464342'
f = open(folder+'reactions_'+pmid+'_curated_model.txt','r')
text=f.read().splitlines()
f.close()

F = []
dict_var = {}
n_var=0
for el in text:
    new_var = el.split('\t')[0].split('=')[1].strip(' ').replace('/','_').replace('-','_')
    rule = '('+el.split('\t')[0].split('=')[0].replace('/','_').replace('-','_').replace('!','NOT ').replace('+',' AND ')+')'
    try:
        index = dict_var[new_var]
        F[index] += ' OR '+rule
    except KeyError:
        dict_var.update({new_var:n_var})
        index = n_var
        F.append(new_var + ' = '+rule)
        n_var+=1
          

g = open('update_rules_models_in_literature_we_randomly_come_across/'+pmid+'_curated_model.txt','w')
for line in F:
    g.write(line+'\n')
g.close()

