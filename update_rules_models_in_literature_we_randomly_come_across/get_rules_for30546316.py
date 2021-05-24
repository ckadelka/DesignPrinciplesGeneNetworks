#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:00:16 2020

@author: tbutrie
"""
folder = 'update_rules_models_in_literature_we_randomly_come_across'
pmid = '25163068'

text = """
lasR = IN1
LasR = lasR
rsaL = IN2
lasl = (IN4 OR LasR+) AND (NOT RsaL)
Lasl = lasl
3oxo = Lasl OR IN3
LasR+ = LasR AND 3oxo
lasREG = LasR+
rhlR = IN5 OR LasR+
RhlR = rhlR
rhlI = IN6 OR LasR+ OR RhlR+ OR PqsR+
RhlI = rhlI
C4-HSL = RhlI OR IN7
RhlR+ = RhlR AND C4-HSL
rhlREG = RhlR+
pqsR = LasR+ AND(NOT RhlR+)
PqsR = pqsR
pqsH = LasR+
pqsA-E = PqsR+
phnAB = PqsR+
PQS = pqsH AND pqsA-E AND phnAB
PqsR+ = PqsR AND PQS
pqsRREG = PqsR+
"""

text = text.replace('=','').replace('[','').replace(']','').replace('"','')

textarray = [el.split(', ') for el in text.split('\n')]    

g = open(folder+pmid+'.txt','w')
for el in textarray:
    if len(el)==2:
        g.write((el[0].replace(' ','')+' = '+el[1].replace(',','')).replace('  ',' ')+'\n')
g.close()

#manually delete constant rules after creation, we don't store them in file

