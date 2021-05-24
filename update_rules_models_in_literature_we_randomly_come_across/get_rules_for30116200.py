#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:00:16 2020

@author: tbutrie
"""
folder = 'update_rules_models_in_literature_we_randomly_come_across'
pmid = '30116200'

text = '''[['Myb = NOT(PU1highORPU1lowORHES1ORCD122high)'],
 ['miR155 = MybORFlt3'],
 ['CEBPα = PU1high AND NOTIkaros AND NOTEbf1 AND NOTFlt3'],
 ['EGR1 = PU1high AND NOTGFI1'],
 ['PU1high = NOT(Bcl11bORmiR155ORPax5ORGFI1ORIrf4ORIkarosORId2)AND (EGR1ORPU1lowORPU1high) AND CEBPα'],
 ['PU1low = NOT(Bcl11bORId2)AND (Flt3OR(PU1lowANDIkaros))'],
 ['Flt3 = NOT(HES1ORBcl11bORCD122highORPax5ORPU1high)AND((IkarosANDPU1low)ORFlt3l)'],
 ['Irf4 = NOT(TbetORPU1high)ANDPU1low'],
 ['Ikaros = Irf4 AND NOTPU1high'],
 ['GFI1 = NOTEGR1AND(IkarosOR(HEBORBcl11b)ANDE2A)'],
 ['Foxo1 = NOT(PDK1ORMAPK) AND (Ebf1ORE2A)'],
 ['E2A = Ebf1ANDNOTHEBANDNOTId2'],
 ['Ebf1 = NOTNotch1 AND PU1low AND ((Ebf1ANDE2AANDPax5) OR Il7)'],
 ['Pax5 = Ebf1 AND NOT(Id2ORHEBORPU1high)'],
 ['Bcl6 = (Pax5OREbf1) AND NOTTbet'],
 ['Notch1 = ((Delta AND PU1low) OR (Notch1 AND Bcl11b AND HEB))ANDNOT(Pax5ORRUNX3ORGzmB)'],
 ['Bcl11b = Notch1ANDNOTPU1high'],
 ['HES1 = Notch1ANDNOTPU1high'],
 ['HEB = Notch1 ANDNOTId2'],
 ['CD122low = NOTBcl11b AND (Flt3ORNotch1)'],
 ['CD122high = NOTBcl11b AND (Tbet OR Eomes OR (ETS1ANDRUNX3))'],
 ['Il15R = Il15 AND CD122low'],
 ['PDK1 = Il15R'],
 ['E4BP4 = PDK1ANDNOTBcl11b'],
 ['MAPK = Il15R'],
 ['ETS1 = MAPKOR(NOTPU1lowANDNOTPU1highANDMAPK)'],
 ['Id2 = NOTBcl11b AND ((NOTEbf1 AND NOTGFI1 AND EGR1) OR E4BP4 OR ETS1)'],
 ['Tox2 = ETS1'],
 ['RUNX3 = NOT(Irf4ORBcl11bORBcl6)AND(ETS1ORTbet)'],
 ['Tbet = NOTIkaros AND NOTBcl6 AND NOTFoxo1 AND ((ETS1 AND Tox2) OR Tbet)'],
 ['Eomes = NOTIrf4 AND (E4BP4 OR (E4BP4 AND RUNX3) OR (Eomes AND RUNX3) OR Eomes)'],
 ['Perforin = Eomes'],
 ['GzmB = TbetOREomes'],
 ['IFN_γ = TbetANDRUNX3'],
 ['']]
'''

text = text.replace('/','_').replace('-','_')

textarray = [el.split('\t') for el in text.split('\n')]    

g = open('update_rules_models_in_literature_we_randomly_come_across/'+pmid+'.txt','w')
for el in textarray:
    g.write(el[0]+' = '+el[2]+'\n')
g.close()

