#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 


folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '31949240'

f = open(folder+pmid+'_to_be_processed.txt','r')
text = f.read()
text = text.replace(')',' ) ').replace('(',' ( ').replace('  ',' ').splitlines()
f.close()

rules = []
var = dict()
for line in text:
    if ' = ' in line:
        rules.append(line)
    elif ' : ' in line:
        var.update({line.split(' : ')[0] : line.split(' : ')[1]})

for i in range(len(rules)):
    rulesplit = rules[i].split(' ')
    for j in range(len(rulesplit)):
        try:
            rulesplit[j] = var[rulesplit[j]]
        except KeyError:
            pass
    rules[i] = ' '.join(rulesplit)


g = open(folder+pmid+'.txt','w')
for rule in rules:
    g.write(rule+'\n')    
g.close()

