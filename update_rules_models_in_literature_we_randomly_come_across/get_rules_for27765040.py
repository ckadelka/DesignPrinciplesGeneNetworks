#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 

text = '''1	PTEN	Src	RhoA	ROCK	PIP3	PTEN		Activated ROCK binds and phosphorylates PTEN, which localizes with RhoA and ROCK to the back and sides of polarized, motile cells.								
		0	0	0	0	0		PTEN is recruited to and activated by the small GTPases only at places containing PIP3 for PTEN to act on.								
		0	0	0	1	0		Hence, the activation condition for PTEN is 'AND(RhoA,ROCK,PIP3)'.								
		0	0	1	0	0										
		0	0	1	1	0										
		0	1	0	0	0										
		0	1	0	1	0										
		0	1	1	0	0										
		0	1	1	1	1										
		1	0	0	0	0										
		1	0	0	1	0										
		1	0	1	0	0										
		1	0	1	1	0										
		1	1	0	0	0										
		1	1	0	1	0										
		1	1	1	0	0										
		1	1	1	1	0										
																
2	Integrins	Src	Ras	Integrins												
		0	0	0												
		0	1	1												
		1	0	0												
		1	1	0												
																
3	EGFR	PTP_EGFR	EGF	EGFR												
		0	0	0												
		0	1	1												
		1	0	0												
		1	1	0												
																
4	PTP_EGFR	Integrins	EGFR	PTP_EGFR		Enhanced EGFR recycling by α5β1 integrin is represented through the inhibition of PTP_EGFR by active integrins.										
		0	0	0												
		0	1	1												
		1	0	0												
		1	1	0												
																
5	VAV2	PI3K	EGFR	VAV2		Treatment of cells with the PI3K inhibitor LY294002 prior to EGF stimulation inhibits VAV2 exchange activity.										
		0	0	0		Hence, both PI3K and EGFR are required for VAV2 activation.										
		0	1	0												
		1	0	0												
		1	1	1												
																
6	Shc	PTEN	FAK	Src	EGFR	Shc										
		0	0	0	0	0										
		0	0	0	1	1										
		0	0	1	0	1										
		0	0	1	1	1										
		0	1	0	0	1										
		0	1	0	1	1										
		0	1	1	0	1										
		0	1	1	1	1										
		1	0	0	0	0										
		1	0	0	1	0										
		1	0	1	0	0										
		1	0	1	1	1										
		1	1	0	0	0										
		1	1	0	1	1										
		1	1	1	0	1										
		1	1	1	1	1										
																
7	Src	CSK	FAK	EGFR	Src		Complete inhibition by CSK is assumed.									
		0	0	0	0											
		0	0	1	1											
		0	1	0	1											
		0	1	1	1											
		1	0	0	0											
		1	0	1	0											
		1	1	0	0											
		1	1	1	0											
																
8	CSK	FAK	Src	CSK		Both FAK and Src are required to bring CSK into close proximity of Src by Paxillin binding.										
		0	0	0												
		0	1	0												
		1	0	0												
		1	1	1												
																
9	PI3K	FAK	Src	EGFR	Ras	PI3K										
		0	0	0	0	0										
		0	0	0	1	1										
		0	0	1	0	1										
		0	0	1	1	1										
		0	1	0	0	1										
		0	1	0	1	1										
		0	1	1	0	1										
		0	1	1	1	1										
		1	0	0	0	1										
		1	0	0	1	1										
		1	0	1	0	1										
		1	0	1	1	1										
		1	1	0	0	1										
		1	1	0	1	1										
		1	1	1	0	1										
		1	1	1	1	1										
																
10	FAK	PTEN	Trio	Integrins	FAK	Src	FAK		The sequence of FAK autophosphorylation followed by Src phosphorylation is needed for the full activation of FAK.							
		0	0	0	0	0	0		Also, Trio or integrins is needed for FAK autophosphorylation, hence forming an OR relation.							
		0	0	0	0	1	0		Active FAK is assumed to be able to maintain its autophosphorylated state.							
		0	0	0	1	0	0		Hence, the activation condition for FAK is 'AND(OR(Trio,Integrins,FAK),Src)'.							
		0	0	0	1	1	1									
		0	0	1	0	0	0									
		0	0	1	0	1	1									
		0	0	1	1	0	0									
		0	0	1	1	1	1									
		0	1	0	0	0	0									
		0	1	0	0	1	1									
		0	1	0	1	0	0									
		0	1	0	1	1	1									
		0	1	1	0	0	0									
		0	1	1	0	1	1									
		0	1	1	1	0	0									
		0	1	1	1	1	1									
		1	0	0	0	0	0									
		1	0	0	0	1	0									
		1	0	0	1	0	0									
		1	0	0	1	1	0									
		1	0	1	0	0	0									
		1	0	1	0	1	0									
		1	0	1	1	0	0									
		1	0	1	1	1	1									
		1	1	0	0	0	0									
		1	1	0	0	1	0									
		1	1	0	1	0	0									
		1	1	0	1	1	1									
		1	1	1	0	0	0									
		1	1	1	0	1	1									
		1	1	1	1	0	0									
		1	1	1	1	1	1									
																
11	Sos	Erk	PIP3	Shc	EGFR	Sos		PIP3 relieves autoinhibition of Sos by PIP2.								
		0	0	0	0	0		Shc:Grb2:Sos complex gets localized to the membrane through the interaction of Shc with the phosphorylated receptor.								
		0	0	0	1	0		Hence, PIP3, Shc and EGFR are all required for Sos activation.								
		0	0	1	0	0										
		0	0	1	1	0										
		0	1	0	0	0										
		0	1	0	1	0										
		0	1	1	0	0										
		0	1	1	1	1										
		1	0	0	0	0										
		1	0	0	1	0										
		1	0	1	0	0										
		1	0	1	1	0										
		1	1	0	0	0										
		1	1	0	1	0										
		1	1	1	0	0										
		1	1	1	1	0										
																
12	Ras	Sos	Ras													
		0	0													
		1	1													
																
13	Raf	Erk	Ras	Src	PAK	Raf		Ras binding promotes conformational changes that relieve Raf-1 autoinhibition. 								
		0	0	0	0	0		Src family induces phosphorylation on Raf Tyr340 or Tyr341, whereas PAKs can phosphorylate S338 in the cytosol in a RAS-independent manner.								
		0	0	0	1	0		Activation of Raf-1 is a complex process in which phosphorylation of Ser338–Tyr341 is a critical step.								
		0	0	1	0	0		Hence, the activation condition for Raf is 'AND(Ras,Src,PAK)'.								
		0	0	1	1	0										
		0	1	0	0	0										
		0	1	0	1	0										
		0	1	1	0	0										
		0	1	1	1	1										
		1	0	0	0	0										
		1	0	0	1	0										
		1	0	1	0	0										
		1	0	1	1	0										
		1	1	0	0	0										
		1	1	0	1	0										
		1	1	1	0	0										
		1	1	1	1	0										
																
14	Mek	Raf	Mek													
		0	0													
		1	1													
																
15	Erk	Mek	Erk													
		0	0													
		1	1													
																
16	PIP3	PTEN	PI3K	PIP3												
		0	0	0												
		0	1	1												
		1	0	0												
		1	1	0												
																
17	RhoGDI	Src	PIP3	RhoGDI												
		0	0	1												
		0	1	0												
		1	0	0												
		1	1	0												
																
18	PAK	PI3K	Cdc42	Rac	PAK											
		0	0	0	0											
		0	0	1	1											
		0	1	0	1											
		0	1	1	1											
		1	0	0	1											
		1	0	1	1											
		1	1	0	1											
		1	1	1	1											
																
19	betaPix	PAK	betaPix													
		0	0													
		1	1													
																
20	Stathmin	Erk	PAK	Stathmin												
		0	0	1												
		0	1	0												
		1	0	0												
		1	1	0												
																
21	GEFH1	FAK	Stathmin	GEFH1												
		0	0	0												
		0	1	1												
		1	0	1												
		1	1	1												
																
22	p190RhoGAP	FAK	Src	p190RhoGAP												
		0	0	0												
		0	1	1												
		1	0	1												
		1	1	1												
																
23	ROCK	RhoA	ROCK													
		0	0													
		1	1													
																
24	ARHGAP22	ROCK	ARHGAP22													
		0	0													
		1	1													
																
25	RhoGAP1	Src	PI3K	RhoGAP1												
		0	0	0												
		0	1	1												
		1	0	1												
		1	1	1												
																
26	Trio	FAK	Trio													
		0	0													
		1	1													
																
27	p190RhoGEF	FAK	p190RhoGEF													
		0	0													
		1	1													
																
28	Fyn	PI3K	Fyn													
		0	0													
		1	1													
																
29	LARG	Fyn	FAK	LARG												
		0	0	0												
		0	1	1												
		1	0	1												
		1	1	1												
																
30	RhoA	RhoGDI	RhoGAP1	p190RhoGAP	LARG	VAV2	p190RhoGEF	Trio	GEFH1	RhoA						
		0	0	0	0	0	0	0	0	0						
		0	0	0	0	0	0	0	1	1						
		0	0	0	0	0	0	1	0	1						
		0	0	0	0	0	0	1	1	1						
		0	0	0	0	0	1	0	0	1						
		0	0	0	0	0	1	0	1	1						
		0	0	0	0	0	1	1	0	1						
		0	0	0	0	0	1	1	1	1						
		0	0	0	0	1	0	0	0	1						
		0	0	0	0	1	0	0	1	1						
		0	0	0	0	1	0	1	0	1						
		0	0	0	0	1	0	1	1	1						
		0	0	0	0	1	1	0	0	1						
		0	0	0	0	1	1	0	1	1						
		0	0	0	0	1	1	1	0	1						
		0	0	0	0	1	1	1	1	1						
		0	0	0	1	0	0	0	0	1						
		0	0	0	1	0	0	0	1	1						
		0	0	0	1	0	0	1	0	1						
		0	0	0	1	0	0	1	1	1						
		0	0	0	1	0	1	0	0	1						
		0	0	0	1	0	1	0	1	1						
		0	0	0	1	0	1	1	0	1						
		0	0	0	1	0	1	1	1	1						
		0	0	0	1	1	0	0	0	1						
		0	0	0	1	1	0	0	1	1						
		0	0	0	1	1	0	1	0	1						
		0	0	0	1	1	0	1	1	1						
		0	0	0	1	1	1	0	0	1						
		0	0	0	1	1	1	0	1	1						
		0	0	0	1	1	1	1	0	1						
		0	0	0	1	1	1	1	1	1						
		0	0	1	0	0	0	0	0	0						
		0	0	1	0	0	0	0	1	0						
		0	0	1	0	0	0	1	0	0						
		0	0	1	0	0	0	1	1	1						
		0	0	1	0	0	1	0	0	0						
		0	0	1	0	0	1	0	1	1						
		0	0	1	0	0	1	1	0	1						
		0	0	1	0	0	1	1	1	1						
		0	0	1	0	1	0	0	0	0						
		0	0	1	0	1	0	0	1	1						
		0	0	1	0	1	0	1	0	1						
		0	0	1	0	1	0	1	1	1						
		0	0	1	0	1	1	0	0	1						
		0	0	1	0	1	1	0	1	1						
		0	0	1	0	1	1	1	0	1						
		0	0	1	0	1	1	1	1	1						
		0	0	1	1	0	0	0	0	0						
		0	0	1	1	0	0	0	1	1						
		0	0	1	1	0	0	1	0	1						
		0	0	1	1	0	0	1	1	1						
		0	0	1	1	0	1	0	0	1						
		0	0	1	1	0	1	0	1	1						
		0	0	1	1	0	1	1	0	1						
		0	0	1	1	0	1	1	1	1						
		0	0	1	1	1	0	0	0	1						
		0	0	1	1	1	0	0	1	1						
		0	0	1	1	1	0	1	0	1						
		0	0	1	1	1	0	1	1	1						
		0	0	1	1	1	1	0	0	1						
		0	0	1	1	1	1	0	1	1						
		0	0	1	1	1	1	1	0	1						
		0	0	1	1	1	1	1	1	1						
		0	1	0	0	0	0	0	0	0						
		0	1	0	0	0	0	0	1	0						
		0	1	0	0	0	0	1	0	0						
		0	1	0	0	0	0	1	1	1						
		0	1	0	0	0	1	0	0	0						
		0	1	0	0	0	1	0	1	1						
		0	1	0	0	0	1	1	0	1						
		0	1	0	0	0	1	1	1	1						
		0	1	0	0	1	0	0	0	0						
		0	1	0	0	1	0	0	1	1						
		0	1	0	0	1	0	1	0	1						
		0	1	0	0	1	0	1	1	1						
		0	1	0	0	1	1	0	0	1						
		0	1	0	0	1	1	0	1	1						
		0	1	0	0	1	1	1	0	1						
		0	1	0	0	1	1	1	1	1						
		0	1	0	1	0	0	0	0	0						
		0	1	0	1	0	0	0	1	1						
		0	1	0	1	0	0	1	0	1						
		0	1	0	1	0	0	1	1	1						
		0	1	0	1	0	1	0	0	1						
		0	1	0	1	0	1	0	1	1						
		0	1	0	1	0	1	1	0	1						
		0	1	0	1	0	1	1	1	1						
		0	1	0	1	1	0	0	0	1						
		0	1	0	1	1	0	0	1	1						
		0	1	0	1	1	0	1	0	1						
		0	1	0	1	1	0	1	1	1						
		0	1	0	1	1	1	0	0	1						
		0	1	0	1	1	1	0	1	1						
		0	1	0	1	1	1	1	0	1						
		0	1	0	1	1	1	1	1	1						
		0	1	1	0	0	0	0	0	0						
		0	1	1	0	0	0	0	1	0						
		0	1	1	0	0	0	1	0	0						
		0	1	1	0	0	0	1	1	0						
		0	1	1	0	0	1	0	0	0						
		0	1	1	0	0	1	0	1	0						
		0	1	1	0	0	1	1	0	0						
		0	1	1	0	0	1	1	1	1						
		0	1	1	0	1	0	0	0	0						
		0	1	1	0	1	0	0	1	0						
		0	1	1	0	1	0	1	0	0						
		0	1	1	0	1	0	1	1	1						
		0	1	1	0	1	1	0	0	0						
		0	1	1	0	1	1	0	1	1						
		0	1	1	0	1	1	1	0	1						
		0	1	1	0	1	1	1	1	1						
		0	1	1	1	0	0	0	0	0						
		0	1	1	1	0	0	0	1	0						
		0	1	1	1	0	0	1	0	0						
		0	1	1	1	0	0	1	1	1						
		0	1	1	1	0	1	0	0	0						
		0	1	1	1	0	1	0	1	1						
		0	1	1	1	0	1	1	0	1						
		0	1	1	1	0	1	1	1	1						
		0	1	1	1	1	0	0	0	0						
		0	1	1	1	1	0	0	1	1						
		0	1	1	1	1	0	1	0	1						
		0	1	1	1	1	0	1	1	1						
		0	1	1	1	1	1	0	0	1						
		0	1	1	1	1	1	0	1	1						
		0	1	1	1	1	1	1	0	1						
		0	1	1	1	1	1	1	1	1						
		1	0	0	0	0	0	0	0	0						
		1	0	0	0	0	0	0	1	0						
		1	0	0	0	0	0	1	0	0						
		1	0	0	0	0	0	1	1	1						
		1	0	0	0	0	1	0	0	0						
		1	0	0	0	0	1	0	1	1						
		1	0	0	0	0	1	1	0	1						
		1	0	0	0	0	1	1	1	1						
		1	0	0	0	1	0	0	0	0						
		1	0	0	0	1	0	0	1	1						
		1	0	0	0	1	0	1	0	1						
		1	0	0	0	1	0	1	1	1						
		1	0	0	0	1	1	0	0	1						
		1	0	0	0	1	1	0	1	1						
		1	0	0	0	1	1	1	0	1						
		1	0	0	0	1	1	1	1	1						
		1	0	0	1	0	0	0	0	0						
		1	0	0	1	0	0	0	1	1						
		1	0	0	1	0	0	1	0	1						
		1	0	0	1	0	0	1	1	1						
		1	0	0	1	0	1	0	0	1						
		1	0	0	1	0	1	0	1	1						
		1	0	0	1	0	1	1	0	1						
		1	0	0	1	0	1	1	1	1						
		1	0	0	1	1	0	0	0	1						
		1	0	0	1	1	0	0	1	1						
		1	0	0	1	1	0	1	0	1						
		1	0	0	1	1	0	1	1	1						
		1	0	0	1	1	1	0	0	1						
		1	0	0	1	1	1	0	1	1						
		1	0	0	1	1	1	1	0	1						
		1	0	0	1	1	1	1	1	1						
		1	0	1	0	0	0	0	0	0						
		1	0	1	0	0	0	0	1	0						
		1	0	1	0	0	0	1	0	0						
		1	0	1	0	0	0	1	1	0						
		1	0	1	0	0	1	0	0	0						
		1	0	1	0	0	1	0	1	0						
		1	0	1	0	0	1	1	0	0						
		1	0	1	0	0	1	1	1	1						
		1	0	1	0	1	0	0	0	0						
		1	0	1	0	1	0	0	1	0						
		1	0	1	0	1	0	1	0	0						
		1	0	1	0	1	0	1	1	1						
		1	0	1	0	1	1	0	0	0						
		1	0	1	0	1	1	0	1	1						
		1	0	1	0	1	1	1	0	1						
		1	0	1	0	1	1	1	1	1						
		1	0	1	1	0	0	0	0	0						
		1	0	1	1	0	0	0	1	0						
		1	0	1	1	0	0	1	0	0						
		1	0	1	1	0	0	1	1	1						
		1	0	1	1	0	1	0	0	0						
		1	0	1	1	0	1	0	1	1						
		1	0	1	1	0	1	1	0	1						
		1	0	1	1	0	1	1	1	1						
		1	0	1	1	1	0	0	0	0						
		1	0	1	1	1	0	0	1	1						
		1	0	1	1	1	0	1	0	1						
		1	0	1	1	1	0	1	1	1						
		1	0	1	1	1	1	0	0	1						
		1	0	1	1	1	1	0	1	1						
		1	0	1	1	1	1	1	0	1						
		1	0	1	1	1	1	1	1	1						
		1	1	0	0	0	0	0	0	0						
		1	1	0	0	0	0	0	1	0						
		1	1	0	0	0	0	1	0	0						
		1	1	0	0	0	0	1	1	0						
		1	1	0	0	0	1	0	0	0						
		1	1	0	0	0	1	0	1	0						
		1	1	0	0	0	1	1	0	0						
		1	1	0	0	0	1	1	1	1						
		1	1	0	0	1	0	0	0	0						
		1	1	0	0	1	0	0	1	0						
		1	1	0	0	1	0	1	0	0						
		1	1	0	0	1	0	1	1	1						
		1	1	0	0	1	1	0	0	0						
		1	1	0	0	1	1	0	1	1						
		1	1	0	0	1	1	1	0	1						
		1	1	0	0	1	1	1	1	1						
		1	1	0	1	0	0	0	0	0						
		1	1	0	1	0	0	0	1	0						
		1	1	0	1	0	0	1	0	0						
		1	1	0	1	0	0	1	1	1						
		1	1	0	1	0	1	0	0	0						
		1	1	0	1	0	1	0	1	1						
		1	1	0	1	0	1	1	0	1						
		1	1	0	1	0	1	1	1	1						
		1	1	0	1	1	0	0	0	0						
		1	1	0	1	1	0	0	1	1						
		1	1	0	1	1	0	1	0	1						
		1	1	0	1	1	0	1	1	1						
		1	1	0	1	1	1	0	0	1						
		1	1	0	1	1	1	0	1	1						
		1	1	0	1	1	1	1	0	1						
		1	1	0	1	1	1	1	1	1						
		1	1	1	0	0	0	0	0	0						
		1	1	1	0	0	0	0	1	0						
		1	1	1	0	0	0	1	0	0						
		1	1	1	0	0	0	1	1	0						
		1	1	1	0	0	1	0	0	0						
		1	1	1	0	0	1	0	1	0						
		1	1	1	0	0	1	1	0	0						
		1	1	1	0	0	1	1	1	0						
		1	1	1	0	1	0	0	0	0						
		1	1	1	0	1	0	0	1	0						
		1	1	1	0	1	0	1	0	0						
		1	1	1	0	1	0	1	1	0						
		1	1	1	0	1	1	0	0	0						
		1	1	1	0	1	1	0	1	0						
		1	1	1	0	1	1	1	0	0						
		1	1	1	0	1	1	1	1	1						
		1	1	1	1	0	0	0	0	0						
		1	1	1	1	0	0	0	1	0						
		1	1	1	1	0	0	1	0	0						
		1	1	1	1	0	0	1	1	0						
		1	1	1	1	0	1	0	0	0						
		1	1	1	1	0	1	0	1	0						
		1	1	1	1	0	1	1	0	0						
		1	1	1	1	0	1	1	1	1						
		1	1	1	1	1	0	0	0	0						
		1	1	1	1	1	0	0	1	0						
		1	1	1	1	1	0	1	0	0						
		1	1	1	1	1	0	1	1	1						
		1	1	1	1	1	1	0	0	0						
		1	1	1	1	1	1	0	1	1						
		1	1	1	1	1	1	1	0	1						
		1	1	1	1	1	1	1	1	1						
																
31	Rac	RhoGDI	RhoGAP1	ARHGAP22	RhoA	Trio	betaPix	VAV2	Rac		Trio can activate Rac ONLY when RhoA is inactive.					
		0	0	0	0	0	0	0	0		RhoA only suppresses the Rac activation by Trio.					
		0	0	0	0	0	0	1	1							
		0	0	0	0	0	1	0	1							
		0	0	0	0	0	1	1	1							
		0	0	0	0	1	0	0	1							
		0	0	0	0	1	0	1	1							
		0	0	0	0	1	1	0	1							
		0	0	0	0	1	1	1	1							
		0	0	0	1	0	0	0	0							
		0	0	0	1	0	0	1	1							
		0	0	0	1	0	1	0	1							
		0	0	0	1	0	1	1	1							
		0	0	0	1	1	0	0	0							
		0	0	0	1	1	0	1	1							
		0	0	0	1	1	1	0	1							
		0	0	0	1	1	1	1	1							
		0	0	1	0	0	0	0	0							
		0	0	1	0	0	0	1	0							
		0	0	1	0	0	1	0	0							
		0	0	1	0	0	1	1	1							
		0	0	1	0	1	0	0	0							
		0	0	1	0	1	0	1	1							
		0	0	1	0	1	1	0	1							
		0	0	1	0	1	1	1	1							
		0	0	1	1	0	0	0	0							
		0	0	1	1	0	0	1	0							
		0	0	1	1	0	1	0	0							
		0	0	1	1	0	1	1	1							
		0	0	1	1	1	0	0	0							
		0	0	1	1	1	0	1	0							
		0	0	1	1	1	1	0	0							
		0	0	1	1	1	1	1	1							
		0	1	0	0	0	0	0	0							
		0	1	0	0	0	0	1	0							
		0	1	0	0	0	1	0	0							
		0	1	0	0	0	1	1	1							
		0	1	0	0	1	0	0	0							
		0	1	0	0	1	0	1	1							
		0	1	0	0	1	1	0	1							
		0	1	0	0	1	1	1	1							
		0	1	0	1	0	0	0	0							
		0	1	0	1	0	0	1	0							
		0	1	0	1	0	1	0	0							
		0	1	0	1	0	1	1	1							
		0	1	0	1	1	0	0	0							
		0	1	0	1	1	0	1	0							
		0	1	0	1	1	1	0	0							
		0	1	0	1	1	1	1	1							
		0	1	1	0	0	0	0	0							
		0	1	1	0	0	0	1	0							
		0	1	1	0	0	1	0	0							
		0	1	1	0	0	1	1	0							
		0	1	1	0	1	0	0	0							
		0	1	1	0	1	0	1	0							
		0	1	1	0	1	1	0	0							
		0	1	1	0	1	1	1	1							
		0	1	1	1	0	0	0	0							
		0	1	1	1	0	0	1	0							
		0	1	1	1	0	1	0	0							
		0	1	1	1	0	1	1	0							
		0	1	1	1	1	0	0	0							
		0	1	1	1	1	0	1	0							
		0	1	1	1	1	1	0	0							
		0	1	1	1	1	1	1	0							
		1	0	0	0	0	0	0	0							
		1	0	0	0	0	0	1	0							
		1	0	0	0	0	1	0	0							
		1	0	0	0	0	1	1	1							
		1	0	0	0	1	0	0	0							
		1	0	0	0	1	0	1	1							
		1	0	0	0	1	1	0	1							
		1	0	0	0	1	1	1	1							
		1	0	0	1	0	0	0	0							
		1	0	0	1	0	0	1	0							
		1	0	0	1	0	1	0	0							
		1	0	0	1	0	1	1	1							
		1	0	0	1	1	0	0	0							
		1	0	0	1	1	0	1	0							
		1	0	0	1	1	1	0	0							
		1	0	0	1	1	1	1	1							
		1	0	1	0	0	0	0	0							
		1	0	1	0	0	0	1	0							
		1	0	1	0	0	1	0	0							
		1	0	1	0	0	1	1	0							
		1	0	1	0	1	0	0	0							
		1	0	1	0	1	0	1	0							
		1	0	1	0	1	1	0	0							
		1	0	1	0	1	1	1	1							
		1	0	1	1	0	0	0	0							
		1	0	1	1	0	0	1	0							
		1	0	1	1	0	1	0	0							
		1	0	1	1	0	1	1	0							
		1	0	1	1	1	0	0	0							
		1	0	1	1	1	0	1	0							
		1	0	1	1	1	1	0	0							
		1	0	1	1	1	1	1	0							
		1	1	0	0	0	0	0	0							
		1	1	0	0	0	0	1	0							
		1	1	0	0	0	1	0	0							
		1	1	0	0	0	1	1	0							
		1	1	0	0	1	0	0	0							
		1	1	0	0	1	0	1	0							
		1	1	0	0	1	1	0	0							
		1	1	0	0	1	1	1	1							
		1	1	0	1	0	0	0	0							
		1	1	0	1	0	0	1	0							
		1	1	0	1	0	1	0	0							
		1	1	0	1	0	1	1	0							
		1	1	0	1	1	0	0	0							
		1	1	0	1	1	0	1	0							
		1	1	0	1	1	1	0	0							
		1	1	0	1	1	1	1	0							
		1	1	1	0	0	0	0	0							
		1	1	1	0	0	0	1	0							
		1	1	1	0	0	1	0	0							
		1	1	1	0	0	1	1	0							
		1	1	1	0	1	0	0	0							
		1	1	1	0	1	0	1	0							
		1	1	1	0	1	1	0	0							
		1	1	1	0	1	1	1	0							
		1	1	1	1	0	0	0	0							
		1	1	1	1	0	0	1	0							
		1	1	1	1	0	1	0	0							
		1	1	1	1	0	1	1	0							
		1	1	1	1	1	0	0	0							
		1	1	1	1	1	0	1	0							
		1	1	1	1	1	1	0	0							
		1	1	1	1	1	1	1	0							
																
32	Cdc42	RhoGDI	RhoGAP1	VAV2	betaPix	Cdc42										
		0	0	0	0	0										
		0	0	0	1	1										
		0	0	1	0	1										
		0	0	1	1	1										
		0	1	0	0	0										
		0	1	0	1	0										
		0	1	1	0	0										
		0	1	1	1	1										
		1	0	0	0	0										
		1	0	0	1	0										
		1	0	1	0	0										
		1	0	1	1	1										
		1	1	0	0	0										
		1	1	0	1	0										
		1	1	1	0	0										
		1	1	1	1	0										'''.splitlines()
        

folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '25780058'

textsplit = [el.split('\t') for el in text]

F = []
dict_var = {}
n_var=0
for i,el in enumerate(textsplit):
    if el[0]==str(n_var+1):
        n_var+=1
        first_empty_index = el.index('')
        n = first_empty_index-3
        rule = []
        for j in range(2**n):
            if textsplit[i+j+1][first_empty_index-1] == '1':
                rule.append(' ( ')
                for k in range(n):
                    if k>0:
                        rule[-1] += ' AND '
                    if textsplit[i+j+1][2+k] == '0':
                        rule[-1] += ' NOT '
                    rule[-1]+= el[2+k]
                rule[-1] += ' ) '
        F.append(el[1]+ ' = '+' OR '.join(rule))
    else:
        continue

          

g = open(folder+pmid+'.txt','w')
for line in F:
    g.write(line+'\n')
g.close()

