#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 


text='''Cytokines (purple nodes)
1. VEGF = NFkappaBs | STATp | NFkappaBp
f1 = x22 | x45 | x50
2. EGF = cJUNp
f2 = x55
3. bFGF = ERKs | ERKp
f3 = x29 | x53
4. TNFalpha = TNFalpha
f4 = x4
5. PDGFBB = SMADs | SMADp
f5 = x19 | x44
6. Thiazolidinedione = Thiazolidinedione
f6 = x6
7. TGFb1 = SMADs | SMADp
f7 = x19 | x44
8. IFNgamma = IFNgamma
f8 = x8
Stellate Cell (Blue)
9. VEGFRs = VEGF
f9 = x1
10. EGFRs  = EGF
f10 = x2
11. FGFRs  = bFGF
f11 = x3
12. TNFRs  = TNFalpha
f12 = x4
13. PDGFBBRs = PDGFBB
f13 = x5
14. PPARgammas = Thiazolidinedione
f14 = x6
15. TGFRs  = TGFb1
f15 = x7
16. IFNGRs = IFNgamma
f16 = x8
17. RASs = ( ( VEGFRs ) | ( EGFRs ) | ( FGFRs ) )
f17 = x9 | x10 | x11
18. PI3Ks = ( ( PDGFBBRs ) | ( RASs ) )
f18 = x13 | x17
19. SMADs = TGFRs
f19 = x15
20. STATs = IFNGRs 
f20 = x16
21. RAFs = RASs
f21 = x17
22. NFkappaBs = ( ( TNFRs ) | ( AKTs ) )
f22 = x12 | x28
23. P38s = ( ( MEKs ) | ( P53s ) )
f23 = x24 | x26
24. MEKs = RAFs
f24 = x21
25. PIP3s = ( ~ PTENs ) & ( PI3Ks )
f25 = (x27+1) & x18
26. P53s = ( ~ MDM2s ) & ( P38s )
f26 = (x32+1) & x23
27. PTENs = P53s
f27 = x26
28. AKTs = PIP3s
f28 = x25
29. ERKs = ( PDGFBBRs ) | ( MEKs )
f29 = x13 | x24
30. AP1s = ERKs
f30 = x29
31. P21s = P53s
f31 = x26
32. MDM2s = ( ( AKTs ) |  ( P53s ) )
f32 = x28 | x26
33. Apos = P53s
f33 = x26
34. Pros = ( ~ P21s ) & ( AP1s )
f34 = (x31+1) & x30
35. Migs = ( AKTs ) & ( AP1s ) & ( Acts )
f35 = x28 & x30 & x36
36. Acts = SMADs & (( ~STATs ) | ( ~PPARgammas )) & ( NFkappaBs | AP1s )
f36 = x19 & ((x20 +1) | (x14+1)) & (x22 | x30)
Pancreatic Cell (cyan)
37. EGFRp = EGF | HER2p
f37 = x2 | x40
38. FGFRp = bFGF
f38 = x3
39. TGFRp = TGFb1
f39 = x7
40. HER2p = HER2p
f40 = x40
41. JAK1p = HER2p
f41 = x40
42. PI3Kp = ( EGFRp ) | ( RASp )
f42 = x37 | x43
43. RASp = ( EGFRp ) | ( FGFRp )
f43 = x37 | x38
44. SMADp = TGFRp
f44 = x39
45. STATp = JAK1p
f45 = x41
46. PIP3p = ( ~ PTENp ) & ( PI3Kp )
f46 = (x52+1) & x42
47. RAFp = RASp
f47 = x43
48. P21p = ( ~ STATp ) & ( ( SMADp ) | ( P53p ) )
f48 = (x45+1) & (x44 | x64)
49. MEKp = RAFp
f49 = x47
50. NFkappaBp = AKTp
f50 = x51
51. AKTp = PIP3p
f51 = x46
52. PTENp = P53p
f52 = x64
53. ERKp = MEKp
f53 = x49
54. E2Fp = ( ~ RBp )
f54 = (x57+1)
55. cJUNp = ( ERKp ) | ( JNKp )
f55 = x53 | x59
56. CyclinDp = ( ~ P21p ) & ( NFkappaBp )
f56 = (x48+1) & x50
57. RBp = ( ~ CyclinDp )
f57 = (x56+1)
58. BCLXLp = ( ~ P53p ) & ( ( NFkappaBp ) | ( STATp ) | ( AKTp ) | ( JNKp ) )
f58 = (x64+1) & (x50 | x45 | x51 | x59)
59. JNKp = MEKp
f59 = x49
60. mTORp = ( ~ cJUNp ) & ( AKTp )
f60 = (x55+1) & x51
61. BAXp = ( ~ BCLXLp )
f61 = (x58+1)
62. Beclin1p = ( ~ BCLXLp ) & ( ~ CASPp )
f62 = (x58+1) & (x66+1)
63. MDM2p = ( ~ E2Fp ) & ( AKTp | P53p )
f63 = (x54+1) & (x51 | x64)
64. P53p = ( ~ MDM2p )
f64 = (x63+1)
65. CyclinEp = ( ~ P21p ) & ( E2Fp ) 
f65 = (x48+1) & x54
66. CASPp = ( ~ NFkappaBp ) & ( ( P53p ) | ( Beclin1p ) | ( BAXp ) )
f66 = (x50+1) & (x64 | x62 | x61)
67. Autp = ( ~ mTORp ) & ( ( NFkappaBp ) | ( Beclin1p ) )
f67 = (x60+1) & (x50 | x62)
68. Apop = CASPp
f58 = x66
69. Prop = ( CyclinEp ) & (( JNKp ) | ( cJUNp ))
f69 = x65 & (x59 | x55)'''

folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '35752283'

rules = []
variables = []
count=0
for line in text.splitlines():
    if '.' not in line:
        continue    
    var,rule = line.replace('\t','').strip().split('=')
    rule = rule.replace('~',' NOT ').replace('&',' AND ').replace('|',' OR ')
    var = var.split('.')[-1].strip()
    rule = rule.lstrip()
    
    print(var,'= ',rule)
    
    if rule==var:
        continue
    
    variables.append(var)
    rules.append(rule)
    

g = open(folder+pmid+'.txt','w')
for var,rule in zip(variables,rules):
    g.write(var+' = '+rule+'\n')    
g.close()

