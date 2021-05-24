#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 


folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '30281473'

text = '''EGFR = EGF; 
EFGR = EGF | HBEGF;
IFGR1AB = IGF;
ERBB2 = NRG1;
JAK1 = (EGFR & ~Lapatinib);
STAT3 = JAK1;
IRS1 = IFGR1AB;
GRB2 = ((EGFR & ~Lapatinib)| (EFGR & ~Lapatinib) | IFGR1AB | (ERBB2 & ~Lapatinib));
RAS = GRB2;
MEKK1 = RAS;
RAF = RAS;
MKK4 = MEKK1;
MEK1 = RAF;
PIK3CA = (RAS | (STAT3 & (~CRY & ~HO3867)) | (IRS1) | (ERBB2 & ~Lapatinib));
JNK1 = MKK4;
ERK1 = MEK1;
PIP3 = ((PIK3CA & ~LY294002) | (~PTEN)); 
PDPK1 = PIP3;
AKT = PIP3;
AMPK = LKB1;
GSK3 = ~(AKT); 
TSC1 = (AMPK | ~(AKT)); 
RHEB = ~TSC1;
mTOR = RHEB;
RPS6KB1 = ((ERK1 & ~CRY) | PDPK1 | (mTOR & ~(Temsirolimus))); 
BAD = ~(RPS6KB1 | (AKT)); 
CCND1 = ~GSK3; 
BCL2 = ((STAT3 & (~CRY & ~HO3867)) & ~BAD);
SRFELK1 = (ERK1 & ~CRY) & RPS6KB1; 
FOSJUN = JNK1 & RPS6KB1;
SRFELK4 = (ERK1 & ~CRY) & RPS6KB1; 
SP1 = (ERK1 & ~CRY);'''.replace(' & ',' AND ').replace('~',' NOT ').replace('|',' OR ').replace(';','')

g = open(folder+pmid+'.txt','w')
g.write(text)    
g.close()

