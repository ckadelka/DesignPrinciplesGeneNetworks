#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 

text = """ABA or PA activates SK 

SK activates S1P  
S1P activates GPA1
GPA1 or NO or CA activates PLD 


PLD activates PA 
PA activates RCN1 

ABA activates PYR  
ABA activates CuAO  
ABA or PA activates PI4P  


ABA activates PI3P
PI4P activates PIP2 

ABA & CA
activates PLC 


PLC & PIP2 
activates InSP3 

  InSP3 &               (no MRP5) activates InSP6


SnRK2 activates pH 
 (no PYR & no PA & no ROS & no ATGPX3 & (pH  or ROP11))
 or
 (no PYR & no PA & no ROS & no ATGPX3 & ROP11) 
or
no PYR
or
no PYR & ABH1
activates PP2C









ERA1 activates ROP11
PP2C deactivates SnRK2
(ROS or (CDPK & CA)) & no PP2C activates MAPK



CA & no PP2C activates CDPK

(SnRK2 or (CDPK & CA)) & RCN1 & PA & (PI3P or PI4P) & pH activates RbOH








((CuAO or RBOH) & no ATGPX3) 
or
(CuAO & RBOH)
activates ROS



ROS activates GPX3

ROS & CaM activates NO

CGMP activates CADPR
NO activates CGMP

((SnRK2 or (CDPK & no-PP2C)) & MAPK & (no-ERA1 or MRP5 or no-ABH1)) & no-MALATE activates SLAC1







((no-NO or ROS or pH) & DEPOLAR) activates GORK




(SnRK2 & DEPOLAR) or MALATE activates QUAC

CA activates ATALMT6
SnRK2 or CDPK deactivates KAT1

(CA & no-pH) activates TPK1

PEPC activates Malate
(PP1 & no CA & no pH) activates AHA1




((Icas or CADPR or Ica or InSP6)  &  no CaATPASE  & no CAX ) activates CA





CaM activates CaATPASE
CA activates CaM
((no ABH1 or no ERA1 or MRP5) & ROS & no DEPOLAR) activates Ica





CA activates CBL
CIPK activates CAX1
ABA or InSP6 activates SCAB1

ACTIN activates Icas
CBL & no
PP2C activates CIPK

PA deactivates PP1
(CA & no PI3P & no PI4P) activates ABPS


((ARP23 & SCAB1) & no AtRAC1 & (ROS or (CA & ABPS))) activates ACTIN  





MALATE inhibits PEPC
PP2C activates ATRAC1
(PIP2 or ROS) activates ARP23


(CA or TPK1 or no AHA1 or SLAC1 or QUAC) activates DEPOLAR





(GORK & ACTIN & no MALATE & SLAC1) activates CLOSURE
""".split('\n')


folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '30518777'

rules=[]
rule=''
for line in text:
    if line=='':
        continue
    if 'deactivates ' in line or 'activates ' in line or 'inhibits ' in line:
        if 'deactivates ' in line:
            front,back = line.split('deactivates ')
            rule+=' '+front
            rule = 'NOT ( '+rule+' )'
        elif 'activates ' in line:
            front,back = line.split('activates ')
            rule+=' '+front 
        elif 'inhibits ' in line:
            front,back = line.split('inhibits ')
            rule+=' '+front
            rule = 'NOT ( '+rule+' )'
        rules.append((back+' = '+rule).replace('&',' AND ').replace('+','p').replace('or',' OR ').replace('no-','NOT ').replace('no','NOT').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('-','_'))
        rule=''    
    else:
        rule+=line
        
g = open(folder+pmid+'.txt','w')
for line in rules:
    g.write(line+'\n')
g.close()

