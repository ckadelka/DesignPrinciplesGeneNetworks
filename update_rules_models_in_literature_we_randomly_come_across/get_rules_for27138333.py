#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 


text='''EGF-EGFR    
EGFR-Shc    
EGFR-Grb2    
Shc-Grb2    
EGFR-Pi3k    
Grb2-Gab1    
EGFR-Nck    
Gab1-Pi3k    
Grb2+~p90rsk-Sos1    
Grb2+~Erk12-Sos1    
Pi3k+Ship2+Pip2+~Pten-Pi34p2    
Pi3k+Pip2+~Ship2+~Pten-Pip3    
EGFR+Pi34p2-Vav2    
EGFR+Pip3-Vav2    
Vav2-Rac1    
Sos1+Pi3k+Pip3+Eps8-Sos1E    
Sos1E-Rac1    
Pi34p2+Pdk1+Mtor+~Pp2a-Akt    
Pip3+Mtor+Pdk1+~Pp2a-Akt    
Akt-RacGap1    
Akt+RacGap1-IQGap1    
IQGap1+RacGap1-pRacGap1    
Vav2-RhoA    
Gab1-Shp2    
Gab1+~Shp2-RasGap    
Sos1+~RasGap-Ras    
Ras+Csrc-Raf1    
Ras+Pak1-Raf1    
Raf1-Mek12    
Mek12-Erk12    
Erk12+Pdk1-P90rsk    
Nck+Rac1-Pak1    
~pRacGap1-Rac1    
Sos1-Abi    
Abi-Eps8    
Sos1-Hras    
Hras-Ralgds    
Ralgds-Rap1a    
Rap1A-Ralb    
Ralb-RalbP1    
RalbP1-Rac1    
Pip3-Gab1    
~Rac1-RhoA    
~RhoA-Rac1    
Rac1-Mekk1    
Mekk1-Mek12'''

folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '27138333'

rules = []
variables = dict()
count=0
for line in text.splitlines():
    rule,var = line.replace('\t','').strip().split('-')
    rule = rule.replace('~',' NOT ').replace('+',' AND ')
    try:
        index = variables[var]
        rules[index]+=' OR ('+rule+')'
    except KeyError:
        index = count
        variables.update({var:count})
        count+=1
        rules.append('('+rule+')')

rulesRac1 = []
for ruleRac1,table in zip(rulesRac1,['tab1d','tab1e','tab1f']):


g = open(folder+pmid+'.txt','w')
for var,rule in zip(variables,rules):
    g.write(var+' = '+rule+'\n')    
g.close()

