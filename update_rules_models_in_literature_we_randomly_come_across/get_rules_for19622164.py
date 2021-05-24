#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 


folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '19622164'

text='''1	COL1A1.out	= TGFB1 OR FOS OR JUN	*
2	COL1A2.out	= (TGFB1 AND JUND) OR JUN	*
3	ETS1.out	= FOS AND NOT (TGFB1 AND SMAD4 AND NOT SMAD7)	*
4	FOS.out	= (TGFB1 OR (TNF AND NFKB1)) AND NOT (JUN AND FOS)	
5	JUN.out	= (TNF AND JUN) OR (TGFB1 AND SMAD4 AND NOT SMAD7)	
6	JUNB.out	= (TGFB1 AND NFKB1 AND SMAD4 AND NOT SMAD7) OR (TNF AND NFKB1)	
7	JUND.out	= TGFB1 OR TNF OR (JUND AND NOT FOS)	
8	MMP1.out	= JUND AND NOT (TGFB1 AND FOS)	*
9	MMP3.out	= (TNF AND ((ETS1 AND NFKB1) OR JUN)) OR FOS	
10	MMP9.out	= TNF AND ETS1 AND NFKB1	
11	MMP13.out	= (TNF AND (JUN OR NFKB1)) OR FOS	
12	NFKB1.out	= TNF OR ETS1 OR NFKB1	*
13	SMAD4.out	= TGFB1	*
14	SMAD7.out	= (TGFB1 AND SMAD4 AND NOT SMAD7) OR FOS OR JUN	*
15	TIMP1.out	= (TGFB1 AND SMAD4 AND NOT SMAD7) OR (TNF AND (JUN OR JUNB OR JUND OR NFKB1)) OR FOS	
16	TGFB1.out	= FOS OR JUND	
17	TNF.out	= TNF AND ((ETS1 AND NFKB1) OR JUN)	'''



textsplit = [el.split('\t') for el in text.splitlines()]
rules = []
for el in textsplit:
    rules.append(el[1].replace('.out',' ')+el[2])
g = open(folder+pmid+'_TGF_beta1.txt','w')
for rule in rules:
    g.write(rule+'\n')    
g.close()

text2='''1	COL1A1.out	= (FOS OR JUN) AND NOT (TNF AND ETS1 AND NFKB1)	*
2	COL1A2.out	= JUN AND NOT (TNF AND ETS1 AND NFKB1)	*
3	ETS1.out	= (TNF AND (ETS1 OR JUN)) OR FOS	*
4	FOS.out	= (TGFB1 OR (TNF AND NFKB1)) AND NOT (JUN AND FOS)	
5	JUN.out	= (TNF AND JUN) OR (TGFB1 AND SMAD4 AND NOT SMAD7)	
6	JUNB.out	= (TGFB1 AND NFKB1 AND SMAD4 AND NOT SMAD7) OR (TNF AND NFKB1)	
7	JUND.out	= TGFB1 OR TNF OR (JUND AND NOT FOS)	
8	MMP1.out	= (TNF AND ETS1 AND NFKB1) OR FOS OR JUND	*
9	MMP3.out	= (TNF AND ((ETS1 AND NFKB1) OR JUN)) OR FOS	
10	MMP9.out	= TNF AND ETS1 AND NFKB1	
11	MMP13.out	= (TNF AND (JUN OR NFKB1)) OR FOS	
12	NFKB1.out	= TNF AND (ETS1 OR NFKB1)	*
13	SMAD4.out	= TRUE	*
14	SMAD7.out	= (FOS OR JUN) AND NOT (TNF AND NFKB1)	*
15	TIMP1.out	= (TGFB1 AND SMAD4 AND NOT SMAD7) OR (TNF AND (JUN OR JUNB OR JUND OR NFKB1)) OR FOS	
16	TGFB1.out	= FOS OR JUND	
17	TNF.out	= TNF AND ((ETS1 AND NFKB1) OR JUN)'''

textsplit = [el.split('\t') for el in text.splitlines()]
rules = []
for el in textsplit:
    rules.append(el[1].replace('.out',' ')+el[2])
g = open(folder+pmid+'_TGF_alpha.txt','w')
for rule in rules:
    g.write(rule+'\n')    
g.close()
