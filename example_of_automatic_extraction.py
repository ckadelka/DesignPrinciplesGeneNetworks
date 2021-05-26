#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:03:02 2019

@author: ckadelka
"""




A = '''1	ADAM17 := ADAM17;
2	EGF := AT2R+1;
3	EGFR := (ADAM17+1) * (EGF+1)+1;
4	Ras := EGFR;
5	Raf := (AKT+1) * ((HSP90+1) * (Ras+1)+1);
6	MEK := (Raf+1) * (IGF-1+1) * (EGFR+1)+1;
7	ERK := (PAR1+1) * (TNFa+1) * (TGFB+1) * (MEK+1)+1;
8	RSK := ERK;
9	ELK := RSK;
10	Cdc37 := Cdc37;
11	HSP90 := Cdc37;
12	HSP27 := HSP90+1;
13	ERBB2 := (HSP90+1) * (EGF+1)+1;
14	ACK1 := HSP90;
15	MED15 := MED15;
16	SKIP := SKIP;
17	TGFB := (AR+1) * (SKIP+1) * (MED15+1)+1;
18	SMAD4 := TGFB;
19	SMAD2 := TGFB;
20	SMAD3 := (TGFB+1) * (MED15+1)+1;
21	FOXH1 := SMAD2*SMAD4+1;
22	MED1 := (ERK+1) * (AKT+1)+1;
23	SIRT1 := FHL2;
24	FOXO1:= (SIRT1+1)*PTEN;
25	CDKN2A := Bmi1+1;
26	UBE2C := MED1;
27	FHL2 := FHL2;
28	ERG := ERG;
29	AT2R := MTUS1;
30	PTEN := PTEN;
31	FOXA1 := MED1;
32	JNK := MAPK;
33	NKX3.1 := CDKN2A+1;
34	CASP3 := (NKX3.1*CDKN1B+1) * (JNK+1)+1;
35	TRAP6 := TRAP6;
36	PAR1 := TRAP6;
37	MAPK := (ERK+1) * (Androgen+1)* (EGF+1) * (TGFB+1) * (TNFa+1)+1;
38	Bmi1 := PTEN * (E2F1+1)+1;
39	IL-6 := (MAPK+1) * (WHSC1*NFKB +1) * (BMP-6+1)+1 ;
40	AR := (FOXH1* ERBB2* EGFR* FOXO1* CDKN2A * DAB2IP*EBP1* (SIAH2*NCOR)+1) * ((AKT+1) * (Androgen+1) * (EP300+1) * (Cdc37*Vav3+1) * (MED1+1) * (FHL2+1) * (CACNA1D+1) * (NFKB+1) * (Vav3+1) * (IL-6+1)+1);
41	TNFa := NKX3.1;
42	CCND := (AR+1) * (ERBB2+1) * (IL-8+1) * (WHSC1*NFKB +1)+1;
43	TMPRSS2 := TMPRSS2;
44	CACNA1D := TMPRSS2*ERG;
45	TNFSF11 := TGFB;
46	MTUS1 := MTUS1;
47	EBP1 := EBP1;
48	STAT3 := IL-6;
49	JAK2 := (IL-6+1) * SOCS2+1;
50	PKIB := PKIB;
51	PKA-C := PKIB;
52	CXCL1 := NFKB;
53	BMP-6 := NFKB;
54	IKKa := AKT;
55	AKT := (EBP1 * INPP4B * Fkbp5 +1)* ((HSP90+1) * (PKA-C+1) * * (CXCL1+1) * * (TNFa+1) * (PI3K+1)+1);
56	c-Myc := AR;
57	mTOR := AKT;
58	Androgen := CACNA1D;
59	MMP-9 := (MAPK+1) * (TGFB+1) * (PDEF+1)+1;
60	RhoGAP := RhoGAP;
61	EpCAM := EpCAM;
62	PI3K := (EGF+1) * (IL-4+1) * PTEN * (TGFB+1) * (EpCAM+1) * (KAT5+1) * (Vav3+1)+1;
63	Wnt := Wnt;
64	NFKB := (IKKa+1)* ((EGR-1+1) * (TNFSF11+1) * (WHSC1+1) * (IL-6+1) * (TNFa+1) * (HOXB13+1) * (Wnt+1) * (RhoGAP+1) * (AKT+1) +1);
65	SKP2 := (c-Myc+1) * (mTOR+1)+1;
66	PSA := (NFKB+1) * CDKN2A+1;
67	Vav3 := Vav3;
68	KAT5 := KAT5;
69	INPP4B := AR * NCOR;
70	Fkbp5 := AR;
71	EGR-1 := IGF-1;
72	BCL-2 := NFKB；
73	PDEF := (KAT5+1) * HOXB13+1;
74	Bcl-XL := (AR+1) * (BCL-2+1)+1;
75	BAK := Bcl-XL+1;
76	IL-4 := IL-4;
77	EP300 := (IL-4+1) * (E2F1+1)+1;
78	DAB2IP := Wnt;
79	WHSC1 := AR;
80	NCOR := CK2+1;
81	SIAH2 := SIAH2;
82	CDKN1B := SKP2+1;
83	ZnT4 := ZnT4;
84	HOXB13 := ZnT4+1;
85	E2F1 := (CDKN1A+1)  * (IKKa+1)+1;
86	CDKN1A := (RAC1 * HOXB13+1) *( (SMAD3+1) * (SMAD4+1)+1);
87	IGF-1 := (DAB2IP+1) * SOCS2+1;
88	IL-8 := (ERK+1) * (EGR-1+1) * (NFKB+1) +1;
89	CK2 := CK2;
90	SOCS2 := Androgen;
91	RAC1 := RAC1;
'''

f = open('26603105.txt','w')
for el in A.split('\n')[:-1]:
    rule = el.split('\t')[1][:-1]
    f.write(rule.replace(':=','=')+'\n')
f.close()
