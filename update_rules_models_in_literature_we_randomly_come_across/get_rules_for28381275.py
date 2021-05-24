#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 

text = """KRAS	EGF & ~ERK
c-Myc	(((~APC & ~p53 & ~BRCA1) & (β-catenin | ERK)) | (MEK | Akt)) & ~ARF
APC	~p53 & ((PTEN | GSK3β) & ~Wnt)
p53	((~APC & BRCA1) | ~MDM2 ) & ATM
Wnt	~p53
ATM	DNAdamage | (~p53 & E2F1)
ARF	DNAdamage | (((~p53 & ~Rb) & (E2F1 | β-catenin)) & (KRAS | c-Myc))
RAF	(KRAS | PAK) & (~ERK & ~Akt)
MEK	~ERK & ((EGF | RAF | PAK) & ~PP2A)
PAK	KRAS & ~Akt
PP2A	~EGF
ERK	(~PP2A & ~IκB) & MEK
BRCA1	(E2F1 & (~SLUG & ~Snail)) & ATM
Akt	KRAS & (~PP2A & ~PTEN)
β-catenin	(~p53 & ((~APC | ~GSK3β) | IκB)) & (~Wnt & Akt)
Rb	(~MDM2 & ~Caspase-3) & ((p21 | ATM) & (~CyclinE & ~CyclinD))
MDM2	((p53 | NF-κB) & ~Rb) & ((~ATM & ~ARF) | ~CyclinE & Akt)
Bcl-2	((~p53 & (~E2F1 | ~c-Myc)) | NF-κB) | (Akt & ~Caspase-3)
E2F1	(~ARF & ~Rb) & ((ATM | MDM2) & ~Akt)
BAX	~Bcl-2 & p53
IκB	NF-κB & ((PTEN | Akt) & ~p53)
NF-κB	((CyclinD | Snail) & (~PTEN & ~c-Myc & ~E2F1)) & ((~IκB | Akt) & (~β-catenin & ~GSK3β))
PTEN	p53 & (~Snail & ~GSK3β)
GSK3β	(~Akt & ~ERK) | ~Wnt
Snail	((NF-κB | ERK) | ~p53) & ATM
SLUG	(β-catenin | (~p53 & ~MDM2 & ~GSK3β)) & ERK
MMP	((ERK & NF-κB) | β-catenin) & ~E-cadherin
E-cadherin	(~Snail & ~SLUG) | (~MMP & PAK)
CyclinE	E2F1 & (~p21 & ~Rb)
CyclinD	(~c-Myc & (ERK | β-catenin | NF-κB)) & ~p21
p21	p53 & (~c-Myc & ~Akt)
Caspase-3	(E2F1 | (BAX & ~Bcl-2)) & (~p21 & ~Akt)""".replace('NF-κB','NFKB').replace('-','_').replace('~',' NOT ').replace('β','beta').replace('κ','k').replace('&',' AND ').replace('|',' OR ').split('\n')

folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '28381275'
        
g = open(folder+pmid+'.txt','w')
for line in text:
    g.write(line.split('\t')[0]+' = '+line.split('\t')[1]+'\n')
g.close()

