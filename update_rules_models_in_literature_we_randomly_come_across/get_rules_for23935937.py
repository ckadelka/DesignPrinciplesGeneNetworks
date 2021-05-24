#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 

text = """SK = (ABA | PA)
S1P= (SK)
GPA1= (S1P)
PLD= (GPA1 | NO|CA)
GTG= (ABA &  ¬ GPA1)
PA= (PLD | (PLC & PIP2))
RCN1= (PA)
CuAO= (ABA)
PI4P= (ABA|PA)
PI3P= (ABA)
PIP2= (PI4P)
PLC= (ABA & CA)
InSP3= (PIP2 & PLC)
InSP6= (InSP3 & ¬ MRP5)
PYR= (ABA)
pH=(SnRK2  & ¬ AHA1)
PP2C= (¬ PYR &¬ PA &¬ ROS &¬ ATGPX3 & (pH | ROP11))| (¬ PYR &¬ PA &¬ ROS &¬  ATGPX3 & ROP11)|¬ PYR| (¬ PYR & ABH1)
ROP11= (ERA1)
SnRK2= ( ¬ PP2C )
RBOH= (SnRK2 | (CDPK & CA)) & RCN1 & PA  & (PI3P|PI4P) & pH
ROS= ((CuAO | RBOH) & ¬ ATGPX3)|(CuAO & RBOH)
ATGPX3=(ROS)
NO= (ROS & CaM)
CADPR= (CGMP)
CGMP= (NO)
Ica= ((¬ ABH1 |¬ ERA1| MRP5) & ROS &  ¬ DEPOLAR )
SLAC1=((SnRK2|(CDPK & ¬ PP2C ))& MAPK & (¬ ERA1 | MRP5|¬ABH1))& ¬MALATE
MAPK= ( ROS|(CDPK & CA) ) & ¬ PP2C 
CDPK= (CA & ¬ PP2C) 
CA= ( (Icas | CADPR | Ica | InSP6)  &  ¬ CaATPASE  & ¬ CAX ) 
CBL= (CA)
GORK= ((¬ NO | ROS| pH) & DEPOLAR)
CaATPASE= (CaM)
CaM= (CA)
QUAC= (SnRK2 & DEPOLAR)  
CAX= (CIPK)
SCAB1= (ABA | InSP6)
Icas= (ACTIN)
CIPK= (CBL &  ¬ PP2C)
TPK1= (CA & ¬ pH)
PP1= (¬ PA)
AHA1= (PP1 & ¬ CA & ¬ pH)
ABPS= (CA  & ¬ PI3P & ¬ PI4P)
ACTIN= ((ARP23 & SCAB1) & ¬ AtRAC1 & (ROS|(CA & ABPS)))   
PEPC= (¬ MALATE & ¬ ABA)
MALATE= (PEPC & ¬ QUAC &  ¬ ATALMT6)
AtRAC1= (PP2C)
ARP23= (PIP2 | ROS)
DEPOLAR= (CA | TPK1 | ¬ AHA1 | SLAC1 | QUAC)
ATALMT6= (CA)
KAT1= (¬ CDPK & ¬SnRK2 & ¬ DEPOLAR)
CLOSURE= (GORK & ACTIN & ¬ MALATE & SLAC1)""".split('\n')


folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '30518777'
        
g = open(folder+pmid+'.txt','w')
for line in text:
    g.write(line.replace('=',' = ').replace('&',' AND ').replace('|',' OR ').replace('¬',' NOT ').replace('  ',' ').replace('  ',' ').replace('  ',' ')+'\n')
g.close()

