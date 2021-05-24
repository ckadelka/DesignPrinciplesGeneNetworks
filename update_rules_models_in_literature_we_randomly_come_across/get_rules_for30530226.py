#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:00:16 2020

@author: ckadelka
"""
folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '30530226'

text = '''GADD45beta	1: TGFbeta = 1	
ATM_ATR	1: Damage = 1	
CHEK1_CHEK2	1: ATM_ATR = 1	p53 accmulation leading to apoptosis
p53	1: NOT (Mdm2 = 1) AND (p38MAPK = 1 OR CHEK1_CHEK2 = 1)	
SMAD2/3	1: TGFbeta = 1	
PTEN	1: SMAD2/3 = 1	
p15INK4b	1: SMAD2/3 = 1	p15INK4b accumulation leading to cell cycle arrest
p14ARF	1: E2F1 = 1 OR p38MAPK = 1	
PI3K_AKT	1: TGFbeta = 1 AND NOT PTEN = 1	
BCL2	1: PI3K_AKT = 1	
Mdm2	1: p53 = 1 AND NOT p14ARF = 1	
p38MAPK	1:GADD45beta = 1 AND ATM_ATR = 1	p38MAPK activation (leading to cycle arrest or apoptosis)
p16INK4a	1: p38MAPK = 1	
p21	1: (p53 = 1 OR SMAD2/3 = 1) AND NOT PI3K_AKT = 1	
pRB	1: NOT CdkCyclin = 1	Dephosphorylated pRB bound to E2F
E2F1	1: NOT (pRB = 1) AND NOT (p53 = 1)	
CdkCyclin	1: (Cdc25ABC = 1) AND NOT (p21 = 1) AND NOT (p16INK4a = 1) AND NOT (p15INK4b = 1)	Cell cycle progression
Cdc25ABC	1:(CHEK1_CHEK2 = 1) AND NOT (p38MAPK = 1)	
Cell_Cycle_Arrest	1: NOT (CdkCyclin = 1) OR NOT (E2F1 = 1)	
Proliferation	1: (CdkCyclin = 1) AND E2F1 = 1	
Apoptosis	1: p53 = 1 AND NOT BCL2 = 1'''

textarray = [el.split('\t') for el in text.replace('\u202f=\u202f1','').replace('1:','').replace('/','_').split('\n')]

F = []
for el in textarray:
    F.append(el[0]+' = '+el[1])

g = open(folder+pmid+'.txt','w')
for f in F:
    g.write(f+'\n')
g.close()

