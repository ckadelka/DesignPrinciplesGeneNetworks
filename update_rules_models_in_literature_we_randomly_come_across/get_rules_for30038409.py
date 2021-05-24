#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:00:16 2020

@author: ckadelka
"""

text = '''Myb ¬(PU1high∨PU1low∨HES1∨CD122high)
miR155 Myb∨Flt3
CEBPα PU1high ∧ ¬Ikaros ∧ ¬Ebf1 ∧ ¬Flt3
EGR1 PU1high ∧ ¬GFI1
PU1high ¬(Bcl11b∨miR155∨Pax5∨GFI1∨Irf4∨Ikaros∨Id2)∧ (EGR1∨PU1low∨PU1high) ∧ CEBPα
PU1low ¬(Bcl11b∨Id2)∧ (Flt3∨(PU1low∧Ikaros))
Flt3 ¬(HES1∨Bcl11b∨CD122high∨Pax5∨PU1high)∧((Ikaros∧PU1low)∨Flt3l)
Irf4 ¬(Tbet∨PU1high)∧PU1low
Ikaros Irf4 ∧ ¬PU1high
GFI1 ¬EGR1∧(Ikaros∨(HEB∨Bcl11b)∧E2A)
Foxo1 ¬(PDK1∨MAPK) ∧ (Ebf1∨E2A)
E2A Ebf1∧¬HEB∧¬Id2
Ebf1 ¬Notch1 ∧ PU1low ∧ ((Ebf1∧E2A∧Pax5) ∨ Il7)
Pax5 Ebf1 ∧ ¬(Id2∨HEB∨PU1high)
Bcl6 (Pax5∨Ebf1) ∧ ¬Tbet
Notch1 ((Delta ∧ PU1low) ∨ (Notch1 ∧ Bcl11b ∧ HEB))∧¬(Pax5∨RUNX3∨GzmB)
Bcl11b Notch1∧¬PU1high
HES1 Notch1∧¬PU1high
HEB Notch1 ∧¬Id2
CD122low ¬Bcl11b ∧ (Flt3∨Notch1)
CD122high ¬Bcl11b ∧ (Tbet ∨ Eomes ∨ (ETS1∧RUNX3))
Il15R Il15 ∧ CD122low
PDK1 Il15R
E4BP4 PDK1∧¬Bcl11b
MAPK Il15R
ETS1 MAPK∨(¬PU1low∧¬PU1high∧MAPK)
Id2 ¬Bcl11b ∧ ((¬Ebf1 ∧ ¬GFI1 ∧ EGR1) ∨ E4BP4 ∨ ETS1)
Tox2 ETS1
RUNX3 ¬(Irf4∨Bcl11b∨Bcl6)∧(ETS1∨Tbet)
Tbet ¬Ikaros ∧ ¬Bcl6 ∧ ¬Foxo1 ∧ ((ETS1 ∧ Tox2) ∨ Tbet)
Eomes ¬Irf4 ∧ (E4BP4 ∨ (E4BP4 ∧ RUNX3) ∨ (Eomes ∧ RUNX3) ∨ Eomes)
Perforin Eomes
GzmB Tbet∨Eomes
IFN-γ Tbet∧RUNX3'''

text = text.replace('!',' NOT ').replace('|',' OR ').replace('&',' AND ')

textarray = [el.split('\t') for el in text.split('\n')]    

g = open(folder+pmid+'.txt','w')
for el in textarray:
    g.write(el[0]+' = '+el[1]+'\n')
g.close()

