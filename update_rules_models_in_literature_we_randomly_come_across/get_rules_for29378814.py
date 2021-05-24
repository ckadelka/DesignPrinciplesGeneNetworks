#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 


folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '29378814'

text='''Activin_A/Nodal	Cytokine	(SignalACT or Oct4) and (not Sox2) and (not Lefty1) and (not Gbx2)
BMP4	Cytokine	SignalBMP or Gbx2 or Tbx3 or Myc
Dnmt3b	Enzyme	(Mycn or Tcf3 or EpiTFs) and not (Cdx2 or Klf4)
EpiTFs	Lineage TF	(SignalBMP or Pitx2 or Dusp6) and (not Cdx2) and not (Klf4 and Sox2)
Esrrb	Pluripotency TF	(Klf4 or Klf2 or Nanog or (SignalWNT and (not Tcf3))) and (not EpiTFs)
Fgfr2	Receptor	((SignalFGF) or Gcnf or Cdx2) and (not Nanog) and (not Oct4)
Gata6	Lineage TF	(Gata6 or SignalERK) and (not Klf2) and (not Nanog) and (not Fgf4)
Gbx2	Pluripotency TF	((SignalWNT and (not Tcf3)) or SignalLIF) and ((Esrrb or Jarid2) and not (Tbx3))
Gcnf	Lineage TF	(Gata6 or Cdx2) and (not EpiTFs)
Jarid2	Pluripotency TF	Klf4 or Oct4
Klf4	Pluripotency TF	(SignalLIF or ((Klf2 or Klf4) and Nanog and Esrrb and (Oct4 and Sox2))) and (not EpiTFs)
Nanog	Pluripotency TF	(Nanog or SignalACT or (Oct4 and Sox2) or Tbx3 or Lrh1 or Klf4) and not (Tcf3 or Gata6)
Oct4	Pluripotency TF	(((Oct4 and Sox2) or Nanog or Klf2 or Klf4) and not (Cdx2 and Oct4) and (not Dnmt3b or Klf2)) or (((Oct4 and Sox2) or Nanog or Lrh1 or Klf2 or Klf4) and (not Gcnf) and (Dnmt3b and (not Klf2)))
Smad6	Signal antagonist	(SignalBMP or Gata6) and (not Oct4)
Smad7	Signal antagonist	(Oct4 or Nanog or Esrrb or Klf4 or Tbx3) and (not Gbx2) and (not Jarid2)
Sox2	Pluripotency TF	Nanog or (Oct4 and Sox2)
Tcf3	Lineage TF	(Nanog or Oct4) and (not SignalWNT)
Cdx2	Lineage TF	(SignalBMP or Cdx2) and not (Cdx2 and Oct4)
Dusp6	TF	SignalERK
Fgf4	Pluripotency TF/Cytokine	Esrrb or Nanog or (SignalWNT and (not Tcf3))
Klf2	Pluripotency TF	((Sox2 and Klf4) or Mycn) and (not Pitx2) and (not Dusp6)
Lefty1	TF	(SignalACT or (SignalWNT and (not Tcf3)) or Mycn or (Oct4 and Sox2)) and (not Jarid2) and (not Fgf4)
Lrh1	Pluripotency TF	(Tbx3 or Klf4 or (Oct4 and Sox2)) and (not Tcf3)
Mycn	TF	(Oct4 and Sox2) and (not Nanog)
Pitx2	TF	(SignalACT or (SignalWNT and (not Tcf3))) and (not Sox2) and (not Jarid2)
Tbx3	Pluripotency TF	(SignalPI3K or Tbx3) and (Esrrb or Nanog or Klf4) and (not SignalERK) and (not Tcf3)
Myc	Pluripotency TF	((SignalERK or (SignalWNT and (not Tcf3))) or SignalLIF or Gbx2) and (not Nanog)
Pecam1	Pluripotency Marker	(Klf2 or Nanog) and (not EpiTFs)
Rex1	Pluripotency Marker	(Nanog or Sox2 or Lrh1 or Klf2 or Esrrb) and (not EpiTFs)'''

text = text.replace(')',' ) ').replace('(',' ( ').replace(' or ',' OR ').replace(' and ',' AND ').replace(' not ',' NOT ').replace('  ',' ').replace('  ',' ').replace('/','_')

textsplit = [el.split('\t') for el in text.splitlines()]
rules = []
for el in textsplit:
    rules.append(el[0]+' = '+el[2])
    
g = open(folder+pmid+'.txt','w')
for rule in rules:
    g.write(rule+'\n')    
g.close()
