#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 

text = """ER	ER	OS	OS
FasL	FasL	TNF α	TNF α or NF κB
IL-6	IL-6 or NF κB		
GRP78	ATF6 or XBP1 or ATF4 and (not ER)	ATF6	GRP78
PERK	GRP78 and (not DNAJC3)	IRE1	BAX or BAK or GRP78
EIF2S1	GADD34 and (not PERK)	DNAJC3	ATF6 or XBP1
ATF4	EIF2S1	CHOP	ATF6 or ATF4
XBP1	IRE1	GADD34	CHOP
TNFR1	TNF α	TNFR2	TNF α
TRAF2	IRE1 or TNFR2 or TRADD	ASK1	OS or TRAF2 or DAXX
JNK	OS or ASK1 or GADD45	p38	OS or ASK1
BCL2	(not JNK) and (not CHOP) and (not P53) and (not BAD)	BID	CASP8 and (not BCL2)
BAX	JNK or P53 and (not BCL2)	BAK	BAX and (not BCL2)
DIABLO	BAX or BAK or BID	HtrA2	BAX or BAK or BID
FasR	FasL	TRADD	TNFR1
DAXX	FasR	RIPK1	FasR or TRADD
RAIDD	RIPK1	FADD	FasR or TRADD
CASP8	RAIDD or FADD or CASP3 or CASP6	CASP9	RAIDD or CASP8 or CASP3 or APAF1 or CASP12 and (not XIAP) and (not AKT)
CASP3	CASP9 or CASP8 and (not XIAP)	CASP7	CASP9 or CASP8 or CASP3 or CASP6 and (not XIAP)
CASP6	CASP7 or CASP3		
XIAP	(not DIABLO) and (not HtrA2) and (not CASP3)	CytochromeC	BAX or BAK or BID
APAF1	CytochromeC or P53	Apoptosis	CASP3 or CASP6 or CASP7
INS	INS	INSR	INS
IRS	INSR and (not SOCS3) and (not JNK) and (not IKK β) and (not S6K)	PI3K	IRS or JAK
PIP3	PI3K	PDK1	PIP3
AKT	PDK1 or mTORC2 and (not TRB3)	AS160	AKT
PKC α	PDK1	GLUT4	AKT or AS160 or PKC α
GSK3 β	not AKT	GS	not GSK3 β
FOXO1	PERK and (not AKT) and (not XBP1)	PGC1 α	FOXO1
PEPCK	FOXO1	G6PC	FOXO1
PPAR α	PGC1 α	TRB3	PPAR α or CHOP
TSC1/2	(not AKT) and (not IKK β)	Rheb	not TSC1/2
mTORC1	Rheb	S6K	mTORC1
mTORC2	not S6K	BAD	JNK and (not AKT)
JAK	IL-6 and (not SOCS3)	STAT3	JAK
SOCS3	STAT3	IKK β	TRAF2
NF κB	not IKB α	IKB α	not IKK β""".replace(' β','_beta').replace(' α','_alpha').replace(' κ','_kappa').replace(')',' ) ').replace('(',' ( ').replace('-','_').replace('/','_').replace('not',' NOT ').replace('and',' AND ').replace('or',' OR ').split('\n')

folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '30953496'
        
g = open(folder+pmid+'.txt','w')
for total_line in text:
    total_line = total_line.split('\t')
    for line in [total_line[:2],total_line[2:]]:
        if line[1].strip() in [line[0].strip(),'1','0']:
            print('detected constant: '+'\t'.join(total_line))
        else:
            g.write(line[0]+' = '+line[1]+'\n')    
g.close()

