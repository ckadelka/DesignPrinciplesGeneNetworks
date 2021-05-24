#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:00:16 2020

@author: ckadelka
"""
folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '30024932'

text = '''ATM	1	IR AND (NOT Wip1 OR E2F1)	
ATR	1	IR	l
miR-449a	1	IR	Activation by IR
Sirt-1	1	E2F1 OR NOT miR-449a	
p53-MAIN	1	(ATM OR ATR OR (E2F1 AND 14-3-3s)) AND NOT Mdm2	p53 activation
Mdm2	1	(NOT Wip1 OR p53-MAIN OR RB) AND NOT ATM AND NOT ATR	
p53-Arrest	1	(p53-MAIN OR NOT p53-INP1) AND NOT p53-Killer AND NOT Sirt-1	p53 phosphorylated at Ser-15 and Ser-20
p53-Killer	1	NOT p53-Arrest AND (NOT Sirt-1 OR NOT Wip1) AND p53-MAIN	p53 phosphorylated at Ser-15, Ser-20 and Ser-46
p53-INP1	1	p53-Arrest OR p53-Killer	Control of p53 accumulation
Wip1	1	p53-Arrest	
p21	1	(p53-Arrest OR p53-Killer) AND NOT c-Myc	
14-3-3s	1	p53-Arrest OR p53-Killer	
c-Myc	1	(E2F1 OR NOT RB) AND NOT miR-449a	
E2F1	1	(NOT RB AND ((ATM AND ATR AND NOT miR-449a) OR NOT Sirt-1)) OR c-Myc OR Cdc25ABC	
RB	1	NOT Mdm2 AND NOT Cdc25ABC AND NOT Sirt-1	Dephosphorylated RB bound to E2F1
Cdc25ABC	1	(NOT miR-449a OR c-Myc) AND NOT ATM AND NOT ATR AND NOT 14-3-3s	
Cdc2-CycB	1	Cdc25ABC OR (NOT p21 AND NOT 14-3-3s)	
Proliferation	1	NOT p53-MAIN AND (Cdc2-CycB OR E2F1)	p53-MAIN inhibition and activation of cycle regulators
G2/M-Arrest	1	p21 OR 14-3-3s	G2/M checkpoint arrest phenotype
G2/M-Apoptosis	1	p53-Killer	Apoptosis phenotype'''

text = text.replace('/','_').replace('-','_')

textarray = [el.split('\t') for el in text.split('\n')]    

g = open('update_rules_models_in_literature_we_randomly_come_across/'+pmid+'.txt','w')
for el in textarray:
    g.write(el[0]+' = '+el[2]+'\n')
g.close()

