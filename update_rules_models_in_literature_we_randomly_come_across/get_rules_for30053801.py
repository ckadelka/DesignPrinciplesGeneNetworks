#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:00:16 2020

@author: ckadelka
"""
folder = 'update_rules_models_in_literature_we_randomly_come_across'
pmid = '30053801'

text = '''CmtR	= (Lsr2 OR Rv0081 OR Rv0676c OR CmtR) AND (NOT Rv0324)
CsoR =	(Rv0081 OR CsoR) AND (NOT Rv2324)
DevR =	DevS OR (PhoP AND Rv2034) OR (DevR AND (DevS OR (PhoP AND Rv2034)))
DevS =	DevR
HrcA =	Rv0324
KmtR =	KmtR AND (NOT Rv2034)
KstR =	(Rv0348 OR Rv2324 OR KstR OR (Lsr2 AND Rv0324)) AND (NOT Rv0081)
Lsr2 =	(Rv1985 OR CmtR OR Rv2034 OR KstR OR Lsr2) AND (NOT(PhoP OR Rv0324 OR Rv0474))
PhoP =	(WhiB4 OR Rv0081 OR (Lsr2 AND Rv2034)) AND (NOT Rv2011)
Rv0081 =	DevR AND (NOT (Rv1990c OR Rv0081))
Rv0324 = (Rv0324 OR (Lsr2 AND Rv1985c)) AND (NOT Rv0081)
Rv0348 = Rv0081 AND (NOT (Rv0324 OR Rv1985c OR Rv2034))
Rv0474 =	(Rv1985c OR Rv0081 OR Rv0474) AND (NOT CmtR)
Rv0494 =	(Rv0081 OR Rv0494) AND (NOT Rv0767c)
Rv0653c =	Lsr2 OR Rv0081 OR Rv0653c
Rv0767c =	external_or_internal_regulator
Rv1255c =	NOT Rv1990c
Rv1985c =	(Rv0081 OR Rv1985c) AND (NOT (DevR OR CmtR OR Rv0324))
Rv1990c =	Lsr2 OR Rv0081 OR Rv1990c
Rv2011c =	Rv2011c
Rv2021c	= Rv0081
Rv2034 =	Rv0081 AND (NOT DevR)
Rv2324 =	external_or_internal_regulator
WhiB4 =	Rv0081 AND (Not Rv2011c)
'''

text = text.replace('/','_').replace('-','_')

textarray = [el.split('\t') for el in text.split('\n')]    

g = open('update_rules_models_in_literature_we_randomly_come_across/'+pmid+'.txt','w')
for el in textarray:
    g.write(el[0]+' = '+el[2]+'\n')
g.close()

