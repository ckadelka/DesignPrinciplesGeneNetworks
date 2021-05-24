#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

import numpy as np
import os 


folder = 'update_rules_models_in_literature_we_randomly_come_across/'
pmid = '25063553'

text = '''1	→ ins [reacs_ins]	1	1.0	logical model input
2	ins = ir [reacs_ir]	1	1.0	insulin-to-receptor binding (α-subunit) induces conformational changes that stimulate derepression and intrinsic IR tyrosine kinase activity (β-subunit) leading to receptor auto-/ transphosphorylation and initiation of downstream phosphorylation cascades (for review see [58, 126]); cPKCβ was shown to down-regulate IR autophosphorylation acivity in insulin +/- TPA-treated HEK293 cells, though molecular mechanisms not yet resolved11; cPKCα and nPKCδ, however, seem to function via IRS1 phosphorylation / regulation60 (see below)
3a	irs1_pst = tdum_irs1_pst_irs1_py [reacs_tdum_irs1_pst_irs1_py]	2	0.8	according to common understanding, insulin induces positive-regulatory IRS1-Y phosphorylation via the IR121, 120 and alternatively IGF1R (as shown for primary murine hepatocytes103); however, IRS2 activation seems to substantially depend on insulin-stimulated IR activity40 (as shown for primary murine hepatocytes) and to predominantly account for metabolic and growth-promoting functions of insulin in (murine) liver103; while basal IRS-S/T phosphorylation (kinases not yet unequivocally identified) seems to be required for optimal (or at least permits) subsequent IRS-Y phosphorylation by insulin or IGF1 receptor kinases, the latter, as well as IRS domain functions, were shown to be negatively modulated by hyper-IRS-S/T-phosphorylation40 (for review see [108, 126]); kinases catalysing IRS2-S/T phosphorylation seem to differ with respect to those regulating IRS1 (pointing to isozyme specificity70) and have not yet been unequivocally resolved; IGFR and IRS2 involvement omitted for model reduction purposes
3	ir · !irs1_pst = irs1_py [reacs_irs1_py]
	1	0.8	
4	irs1_py = erk12 [reacs_erk12]	1	0.8	Grb2:SOS (interacting with phospho-IRS and Shc and activating Ras in response to insulin stimulation113) initiates the insulin-induced MAPK pathway, resulting in MEK1/2 and subsequent ERK1/2 activation (for review see [126]; MEK-dependent ERK1/2-T202/Y204 phosphorylation demonstrated for insulin-treated T47D cells3; intermediate effectors not yet integrated for reasons of model reduction); Grb2 (constitutively linked to SOS) specifically interacts with IRS1-pY895 or IRS2-pY911 in response to insulin (as shown for 32-D cells84 and murine brown adipocytes131 or in murine liver127, respectively; IRS1/2:Grb2 interaction also demonstrated using murine muscle cell line C2C1242; for review see [85]), likely initiating MAPK cascade activation; explicit Grb2:SOS involvement omitted / bypassed for model reduction purposes
5	ir · npkcd_g = npkcd [reacs_npkcd]	1	0.8	DAG, generated via PLD-catalysed phosphatidylcholine (PC) hydrolysis or de novo PA / DAG synthesis involving GPAT, allosterically activates nPKCδ via membrane recruitment (for review see [35, 117]); additional Y-phosphorylation of the catalytic domain is supposed to be a specific regulatory mechanism for nPKCδ and to increase its enzymatic activity (reviewed in [117]); supportingly, insulin stimulation was shown to IR-dependently induce Src:nPKCδ interaction and Src-mediated nPKCδ-Y phosphorylation / activation within 1-15 min (as demonstrated for murine skeletal muscle cells101; nPKCδ-Y187/311/332 phosphorylation also shown for insulin-treated murine hepatocytes14); explicit DAG and Src involvement omitted / bypassed for model reduction purposes
6	ir · irs1_py = pdk1 [reacs_pdk1]	1	0.8	insulin-stimulated IR directly interacts with PDK1 (therefore anchoring membrane-targeted PDK1 at the cell surface) and phosphorylates PDK1-Y9, likely depending on IR C-terminal Y1354/1360 as well as prior PI3K / PIP3-mediated PDK1 membrane recruitment and leading to subsequent PDK1-Y373/376 (auto)phosphorylation (Hsp90 / Src involvement as shown for pervanadate-treated HEK293 cells141 cannot be excluded) and activation (as 
demonstrated for insulin-treated L6 skeletal-muscle cells37); PIP3-regulated plasma membrane translocation might ensure appropriate conformational changes and localisation necessary for optimal PDK1 activation and activity (for review see [129]); IRS1/2 associate with regulatory PI3K subunit p85 (likely due to IRS-pY:p85-SH2 interaction42, 121; relevance of IRS1-Y612 and Y632 demonstrated for insulin-treated NIH-3T3 fibroblasts33) in response to insulin (as shown for primary murine hepatocytes103), relieving allosteric PI3K inhibition and inducing catalytic p110 activity (for review see [126]); explicit PI3K and PIP3 involvement omitted / bypassed for model reduction purposes
7a	pp2a = tdum_pp2a_akt [reacs_tdum_pp2a_akt]	2	0.4	insulin induces positive-regulatory Akt-T308/S473 phosphorylation (as shown for HepG2 cells9) generally requiring prior PIP3-mediated Akt and PDK1 membrane recruitment via PH domains; PDK1-catalysed Akt-T308 (kinase domain / T-loop) phosphorylation (as shown for insulin-treated CHO cells137) subsequently initialises Akt activation (for review see [133]); mTORC2 has been shown to phosphorylate Akt-S473 (C-terminal regulatory domain) in various human cancer cell lines (as demonstrated by knockdown studies106), causing full Akt activation (for review see [133]); however, relevance for hepatic insulin signalling has to be checked; TRB3 binds to and interferes with insulin-induced Akt-S473 phosphorylation (as shown for insulin-treated HepG2 cells31); PP2A seems to mediate Akt-S473/T308 dephosphorylation within the context of hyperosmotic stress (as shown for sorbitol-treated HEK293 cells78); but: influence so far regarded as secondary as relevance for hepatic insulin signalling has not yet been confirmed; mTORC2 involvement omitted for model reduction purposes
7b	trb3 = tdum_trb3_akt [reacs_tdum_trb3_akt]	2	1.0	
7	pdk1 ∙ !pp2a · !trb3 = akt [reacs_akt]	1	0.6	
8	!akt · !pka · basal_act = gsk3 [reacs_gsk3]	1	0.6	activated Akt likely catalyses the inhibitory GSK3α-S21 / GSK3β-S9 phosphorylation (as detected in insulin-treated murine hepatocytes86) in response to insulin (causal relation proved for L6 myotubes24); PKA has also been shown to induce inhibitory GSK3α-S21 / GSK3β-S9 phosphor-ylation (as demonstrated by overexpression studies in HEK293 cells), possibly linking GSK3 activity to changes in intracellular cAMP levels34

9	ir · pdk1 = cpkca [reacs_cpkca]	1	0.8	according to traditional PKC activation models, intracellular calcium-to-cPKC-C2 domain binding initially increases cPKC membrane affinity, the latter in turn facilitating DAG:cPKC-C1 domain interaction; subsequent conformational changes (owing to C1 / C2 domain association with membranes), leading to autoinhibitory pseudosubstrate domain displacement, release the cPKC substrate-binding pocket, finally allowing for cPKC-mediated phosphorylation of membrane substrates; additional, recently described, sequential ‘priming’ S/T phosphorylations within the cPKC-activation loop (likely attributable to PDK1 action; complexes with the C-terminus of membrane-localised unphosphorylated PKC) induce cPKC autophosphorylation and are thus essential to generate catalytically competent cPKC (for review see [117]); supportingly, insulin was shown to PI3K-dependently induce cPKCα-to-membrane translocation in rat adipocytes116; as generally accepted, PLCγ calcium-dependently hydrolyses PIP2, thereby generating IP3 and DAG74, the latter in turn likely stimulating PKCα and subsequent PLD activation in response to insulin (based on studies in HEK293 cells114); DAG principally initialises cPKC / facilitates nPKC membrane association (for review see [117]); IP3 induces Ca2+ release from the endoplasmatic reticulum (ER) via IP3 receptors (for review see [46]), resulting in an intracellular Ca2+ increase, inter alia up-regulating C2 domain-related cPKC membrane affinity; PLC- and IP3-dependent mobilisation of intracellular Ca2+ stores demonstrated for insulin-treated primary rat hepatocytes100; insulin (PI3K- and SH2 domain in- but PLCγ-PH-EF domain dependently; as shown for HEK293 and CHO-IR cells69) induces PLCγ-to-IR recruitment and positive-regulatory PLCγ-Y783 phosphorylation (as demonstrated for HepG2 cells and primary rat hepatocytes69); explicit PLCγ, IP3, calcium, and DAG involvement omitted / bypassed for model reduction purposes
10	cpkca = irs1_pst [reacs_irs1_pst]	1	0.8	likely by catalysing negative-regulatory S/T phosphorylation, cPKCα inhibits IRS1 tyrosine phosphorylation (as shown for insulin +/- TPA-treated HEK293 cells60; for review see [108])

11	erk12 = irs1_pst [reacs_irs1_pst]	1	0.4	ERK1/2 mediate negative-regulatory IRS1-S612 phosphorylation (as shown for glucose-treated RIN β-cells6), potentially functioning as feedback inhibitors during insulin signalling; but: relevance for hepatic hormonal signalling has to be checked!
12	npkcd = irs1_pst [reacs_irs1_pst]	1	0.8	nPKCδ directly phosphorylates IRS1 on at least 18 sites (e.g. S307, S323, S574), thereby inhibiting positive-regulatory, IR-mediated IRS1-Y phosphorylation (as demonstrated in vitro and for insulin +/- PMA-treated H4IIE cells41)

13	→ x5p [reacs_x5p]	1	1.0	logical model input / metabolic output 
14	x5p = pp2a [reacs_pp2a]	1	0.8	X5P activates cytosolic as well as nuclear PP2A (as shown in livers of high-carbohydrate-fed rats57); but: the concrete mechanism has yet to be resolved!
15	→ gcg [reacs_gcg]	1	1.0	logical model input
16	gcg = gcgr [reacs_gcgr]	1	1.0	glucagon-to-GCGR binding results in activation of the coupled G proteins, subsequent induction of adenylat cyclases and increasing intracellular cAMP levels, the latter allosterically triggering PKA activation in hepatocytes (for review see [65])

17	gcgr · !akt = pka [reacs_pka]	1	0.8	glucagon-to-GCGR binding results in activation of the coupled G proteins, subsequent induction of adenylat cyclases (for review see [122]) and increasing intracellular cAMP levels, the latter allosterically triggering PKA activation in hepatocytes (for review see [65]); PDE3B (Akt-dependently → insulin-induced Akt catalyses positive-regulatory PDE3B-S273 phosphorylation, leading to PDE3 activation (as shown for insulin-treated 3T3-L1 adipocytes [64])) hydrolyses cAMP, thereby lowering cytosolic cAMP levels (for review see [77]); explicit AC, PDE3B and cAMP involvement omitted / bypassed for model reduction purposes
Effects on transcription (co-)factor activities
18	cpkca = sp1 [reacs_sp1]	1	0.8	insulin rapidly induces SP1 nuclear translocation, nuclear cPKCα:SP1 association and cPKCα-mediated positive-regulatory SP1 serine phosphorylation / activation within 5 min (as shown for L6 skeletal muscle cells48, 47)

19	!pdk1 = cbp [reacs_cbp]	1	0.8	insulin aPKCι / λ-dependently stimulates CBP-S436 phosphorylation, triggering the dissociation of the CREB-CBP-TORC2 transcription complex (as detected in vitro43); insulin-induced (and PI3K-dependent) positive-regulatory aPKCι/λ-T403/555 phosphorylation (as shown for primary mouse hepatocytes43) likely requires PIP3 synthesis and PDK1 activation: the former might interact with residues within or near the aPKC pseudosubstrate site (regulatory domain), causing conformational changes that release the substrate-binding site, increasing PDK1 access to the aPKC activation loop (e.g. aPKCζ-T410; activation loop phosphorylation generally essential to generate catalytically competent c/aPKCs (as reviewed in [117])), and inducing subsequent auto(trans)-phosphorylation (e.g. at aPKCζ-T560)104; PDK1 was shown to up-regulate positive-regulatory aPKCζ-T410 phosphorylation in L6 myotubes in response to insulin104; but: relevance with respect to hepatic aPKCι/λ activation has to be checked; explicit aPKCι and PIP3 involvement omitted / bypassed for model reduction purposes
20a	sik1 = tdum_sik1_torc2 [reacs_tdum_sik1_torc2]	2	0.8	glucagon induces TORC2-S171 dephosphorylation, the latter triggering nuclear TORC2 translocation and co-activator activity (as shown in murine liver68); but: enzymatic effectors not yet clarified; SIK1 counteracts this effect, (re )phosphorylating TORC2-S171 (as shown for glucagon-treated primary murine hepatocytes68)

20	gcgr · !sik1 = torc2 [reacs_torc2]	1	0.8	
21a	cbp = dum_cbp_or_torc2_pgc1a_g [reacs_dum_cbp_or_torc2_pgc1a_g]	1	1.0	phosphorylated / activated CREB and dephosphorylated TORC2 (likely acting alternatively to CREB co-activator CBP68) seem to cooperatively contribute to PGC1α mRNA expression upon glucagon stimulation (as shown for FSK-treated HepG2 cells45 or primary rat hepatocytes68)

21b	torc2 = dum_cbp_or_torc2_pgc1a_g [reacs_dum_cbp_or_torc2_pgc1a_g]	1	0.4	
21	creb · dum_cbp_or_torc2_pgc1a_g = pgc1a_g [reacs_pgc1a_g]	1	0.8	
22	!akt · pgc1a_g = pgc1a [reacs_pgc1a]	1	0.8	Akt catalyses negative-regulatory PGC1α-S570 phosphorylation (as shown for insulin-treated H4IIe rat hepatoma cells72) likely inhibiting PGC1α promoter recruitment
23	ir = creb [reacs_creb]	1	0.8	comparably to glucagon, insulin induces positive-regulatory CREB-S133 phosphorylation in murine liver68; but: enzymatic effectors not yet clarified!
24	pka = creb [reacs_creb]	1	0.4	mediating cAMP effects on gene expression, PKA was shown to catalyse positive-regulatory CREB-S133 phosphorylation39, thereby facilitating CREB interaction with its transcriptional coactivator CBP92 (for review see [111]); but: hepatocellular relevance has to be checked!
25	creb · torc2 = sik1 [reacs_sik1]	1	0.8	glucagon seems to up-regulate SIK1 mRNA levels in a CREB- and TORC2-dependent manner (as shown for forskolin-treated primary murine hepatocytes68), thus initiating TORC2 re-phosphorylation / inactivation within a negative feedback loop
26	!pka · basal_act = hnf4a [reacs_hnf4a]	1	0.6	PKA was shown to phosphorylate HNF4α (possibly at S133/134; as shown in vitro), likely mediating the negative-regulatory effect of glucagon or PKA inducer IBMX on HNF4α DNA-binding activity (as shown in rat liver135)

27	basal_exp = lxra_g [reacs_lxra_g]	1	1.0	fatty acids (including PUFAs) have been shown to up-regulate basal LXRα mRNA and protein levels, likely inducing PPARα-to-PPRE binding within the lxra promoter (as shown in livers of TTA-, Wy14,643-, soy oil- or linolenic acid-treated rats128)

28	ppara = lxra_g [reacs_lxra_g]	2	0.8	
29	→ gluc_in [reacs_gluc_in]	1	1.0	logical model input / metabolic output
30	→ g6p [reacs_g6p]	1	1.0	logical model input / metabolic output
31a	gluc_in = dum_gluc_in_or_g6p_lxra [reacs_dum_gluc_in_or_g6p_lxra]	1	1.0	in addition to (and possibly working in combination with) oxysterols26, 71, D glucose and D-glucose-6-phosphate act as direct LXR agonists, inducing co-activator (SRC1, steroid receptor coactivator-1) recruitment and expression of LXR target genes (as shown for HepG2 cells and in murine liver81); oxysterol involvement omitted for model reduction purposes
31b	g6p = dum_gluc_in_or_g6p_lxra [reacs_dum_gluc_in_or_g6p_lxra]	1	1.0	
31	lxra_g · dum_gluc_in_or_g6p_lxra = lxra [reacs_lxra]	1	1.0	
32a	lxra = dum_lxra_or_basal_exp_insig2a [reacs_dum_lxra_or_basal_exp_insig2a]	1	0.8	in contrast to Insig1, Insig2 was shown to be lowly but constitutively and SREBP-independently expressed (as shown for CHO cells138); insulin specifically stimulates the turnover / degradation of Insig2a mRNA via its 3’-UTR region (as shown for insulin-treated rat hepatocytes143); opposing the effect of insulin, activated LXRα up-regulates Insig2a mRNA expression (as shown for TO-901317-treated rat hepatocytes44) and likely primes the SREBP1c system for insulin-mediated SREBP1c cleavage → ensures lipogenesis only in case of glucose abundance44; but: putative molecular mediators still have to be identified!
32b	basal_exp = dum_lxra_or_basal_exp_insig2a [reacs_dum_lxra_or_basal_exp_insig2a]	1	0.4	
32	!ir · dum_lxra_or_basal_exp_insig2a = insig2a [reacs_insig2a]	1	0.6	
33	basal_exp = pchrebp [reacs_pchrebp]	1	1.0	LXR:RXR heterodimers specifically associate with chrebp promoter elements and up-regulate (constitutive) hepatic ChREBP mRNA expression (as shown in murine liver in response to LXR and / or RXR ligands17)

34	lxra = pchrebp [reacs_pchrebp]	1	0.4	
				
35	pp2a ∙ !pka · pchrebp = nchrebp [reacs_nchrebp]	1	0.8	PKA catalyses negative-regulatory S196/T666 phosphorylation of cytosolic but also nuclear ChREBP (as shown in vitro59), thereby impeding ChREBP nuclear translocation (as demonstrated for glucose-treated primary rat hepatocytes59; possibly owing to phospho-dependent 14-3-3 binding to the ChREBP-NLS motif and CRM1 association, blocking importin α/β interaction and preventing nuclear import / inducing nuclear export105) and DNA binding (as shown in cAMP-treated rat liver59; likely due to impaired ChREBP:Mlx interaction105) in response to glucagon / cAMP; PP2A in turn dephosphorylates PKA target sides, enabling nuclear localisation, subsequent DNA binding, and transactivation (as shown in vitro59; for review see [27])

36a	lxra * sp1 = dum_lxra_sp1_psrebp1c [reacs_dum_lxra_sp1_psrebp1c]	1	0.6	LXRα, NF-Y, Sp1, and SREBP1 somehow contribute to full transcriptional activation of the srebp1c promoter in response to insulin (as shown for primary rat hepatocytes15, 19, 44; for review see [97]); analysis of the mouse srebp1c promoter revealed nSREBP1c to bind to SRE3 in the SRE complex self-regulating its expression (→ positive feedback-loop / autoloop regulation) and NF-Y to be indispensable for basal promoter activity and SREBP-mediated SREBP1c induction4; LXR:RXR was shown to activate the srebp1c promoter and induce precursor expression (as shown for HepG2 cells144; though not seeming sufficient for TF maturation44), thereby indirectly up-regulating downstream lipogenic genes such as  FAS, SCD1, and ACC110; studies in HEK293 cells / skeletal muscle revealed LXREs28 and SP14 to be important for basal promoter activity though not involved in the stimulatory effect of insulin; hence, individual relevances have to be checked; by reducing nuclear LXR:PXR heterodimer formation (possibly owing to RXR competition), PPARα inhibits LXRα:RXR binding to LXREs within the srebp1c promoter, thereby suppressing srebp1c promoter activity (as shown for T0901317 +/- Wy14,643-treated primary rat hepatocytes and HEK293 cells145); NF Y involvement omitted for model reduction purposes
36b	ppara = tdum_ppara_psrebp1c [reacs_tdum_ppara_psrebp1c]	2	0.8	
36	dum_lxra_sp1_psrebp1c · !ppara = psrebp1c [reacs_psrebp1c]	1	0.6	
37	nsrebp1c = psrebp1c [reacs_psrebp1c]	2	0.4	
				
38	!gsk3 · !insig2a · psrebp1c = nsrebp1c [reacs_nsrebp1c]	1	1.0	Insig2a selectively retains SCAP:pSREBP1c complexes in the ER membrane, therefore preventing their ER-to-Golgi transit via COPII vesicles143, proteolytic pSREBP1c intramembrane processing by Golgi-resident proteases, and liberation of transcriptionally active N-terminal fragments (for review see [30]); insulin blocks C-terminal nSREBP1c-S410 phosphorylation by GSK3, the latter initiating subsequent T426 and S430 phosphorylation within its phosphodegron, leading to SCFFbw7 ubiquitin ligase binding and proteasomal degradation of mature nSREBP1c (as shown for HepG2 cells10)

39a	hnf4a = dum_hnf4a_or_ppara_ppara_g [reacs_dum_hnf4a_or_ppara_ppara_g]	1	0.6	HNF4α induces PPARα expression via a DR-1 response element (αHNF4-RE) within the human ppara promoter (as shown by co-transfection studies in HepG2 cells115); by binding to αHNF4-RE (identified as functional PPRE) PPARα was shown to positively autoregulate its own expression (as shown for fenofibric acid-treated primary human hepatocytes115); PGC1α co-activates the ppara gene promoter, an effect that is augmented by LPIN1, the latter interacting with PGC1α (as shown in murine liver36) as well as with the ppara promoter and co-activating PPARα and HNF4α (as shown by co-transfection studies in HepG2 cells36)

39b	ppara = dum_hnf4a_or_ppara_ppara_g [reacs_dum_hnf4a_or_ppara_ppara_g]	2	0.6	
39	dum_hnf4a_or_ppara_ppara_g · pgc1a · lpin1 = ppara_g [reacs_ppara_g]	1	0.6	
40	erk12 · ppara_g = ppara [reacs_ppara]	1	0.6	insulin induces PPARα-S12/21 (located in ligand-independent AF-1 domain) phosphorylation in a likely ERK1/2-dependent manner, thereby (probably owing to co-repressor dissociation) increasing PPARα transcriptional activity (as shown for insulin-treated HepG2 cells56); however, intermediate effectors have to be resolved!
Effects on (metabolic) target gene expression
41	sp1 = npkcd_g [reacs_npkcd_g]	1	0.8	insulin induces SP1 activation, SP1:pkcd promoter interaction and SP1-dependent up-regulation of nPKCδ mRNA transcription and protein levels (as shown for insulin-treated L6 skeletal muscle cells48); but: involvement of potential co-factors (e.g. NF-κB) has to be checked!
42	ppara · pgc1a = trb3 [reacs_trb3]	1	0.6	PPARα targets the trb3 promoter and increases TRB3 gene expression PGC1-dependently; insulin co-treatment reverses the effect (as shown for GW7,647+insulin-treated primary murine hepatocytes and by overexpression studies in primary rat hepatocytes67)

43a	cbp = dum_cbp_or_torc2_g6pc [reacs_dum_cbp_or_torc2_g6pc]	1	0.4	CREB and its co-activators TORC2 (as shown in murine liver upon glucagon stimulation68) and CBP (as shown for FSK-treated primary rat hepatocytes68) bind to the g6pc promoter with CBP and TORC2 seeming to alternatively facilitate G6PC expression via a CREB-TORC or CREB-CBP pathway (two distinct functional CREs within the human g6pc promoter detected109)

43b	torc2 = dum_cbp_or_torc2_g6pc [reacs_dum_cbp_or_torc2_g6pc]	1	0.8	
43	creb · dum_cbp_or_torc2_g6pc = g6pc [reacs_g6pc]	1	0.8	
44	→ f26p2 [reacs_f26p2]	1	1.0	logical model input / metabolic output
45	nchrebp · f26p2 = g6pc [reacs_g6pc]	1	0.8	ChREBP was shown to be F-2,6-P2-dependently (molecular mechanisms not yet resolved) recruited to the g6pc promoter (glucose-induced ChREBP:g6pc-ChoRE binding demonstrated93), up-regulating G6PC gene expression in response to glucose stimulation (as demonstrated for primary rat hepatocytes7)

46	g6pc →	1	0.8	logical model output / metabolic input
47	hnf4a = gys2_g [reacs_gys2_g]	1	0.6	HNF4α (likely contributing to basal GYS2 expression) and PPARα seem to individually transactivate hepatic GYS2 gene expression via two distinct PPREs (PPARα: DR-1int vs. HNF4α: DR-1prom) though negative cross-talk effects between both TFs (possibly owing to competition for common co-activator proteins) referring to maximal transactivation have been observed (as shown by overexpression studies in HepG2 cells76)

48	ppara = gys2_g [reacs_gys2_g]	1	0.6	
49	nsrebp1c * hnf4a = lgk [reacs_lgk]	1	0.8	nSREBP1c mediates insulin-induced liver-specific LGK expression by binding to SREs within the gk promoter (as shown for insulin-treated primary rat hepatocytes63); HNF4α was shown to somehow contribute to insulin-induced gk promoter activity in primary rat hepatocytes possibly by forming a regulatory transcriptional complex with HIF1 (hypoxia-inducible factor 1) and coactivator p300102; for review see [52]

50	lgk →	1	0.8	logical model output / metabolic input
51	nchrebp · hnf4a · cbp = lpk_g [reacs_lpk_g]	1	0.8	nChREBP mediates the stimulatory effect of glucose on LPK gene expression via ChoRE binding (as shown for primary rat and murine hepatocytes54, 140); HNF4α essentially contributes to optimal LPK expression (as shown for glucose +/- insulin-treated primary rat hepatocytes1), binding to the lpk promoter (DR-1 element; as shown in liver of carbohydrate-fed rats1) and likely directing and / or stabilising nChREBP promoter binding via physical interaction1; co-activator CBP is additionally required for glucose-mediated lpk transactivation (for review see [95])

52	pgc1a = pdhk2_g [reacs_pdhk2_g]	1	0.4	PGC1α induces PDHK2 gene expression (as shown by overexpression studies in primary rat hepatocytes73; no apparent PPARα or GCR involvement51); inhibitory effect of insulin on PGC1 gene expression25 might be reflected by the higher basal level and insulin-sensitivity of PDHK2 expression as compared to PDHK451

53a	cbp = dum_cbp_or_torc2_pepck_c [reacs_dum_cbp_or_torc2_pepck_c]	1	0.4	CREB and its co-activators TORC2 (as shown in murine liver upon glucagon stimulation68) and CBP (as shown for FSK-treated primary rat hepatocytes68) bind to the pepck promoter with CBP and TORC2 seeming to alternatively facilitate PEPCK expression via a CREB-TORC or CREB-CBP pathway; nSREBP1c and nSREBP2 were shown to competitively inhibit PEPCK-C gene expression through interaction with SREs within the pepck promoter (as shown by co-transfection studies in HepG2 cells18), but: relevance for hepatic insulin / glucagon signalling has to be checked; for review see [23]; SREBP2 involvement omitted for model reduction purposes
53b	torc2 = dum_cbp_or_torc2_pepck_c [reacs_dum_cbp_or_torc2_pepck_c]	1	0.8	
53	creb · dum_cbp_or_torc2_pepck_c · !nsrebp1c = pepck_c [reacs_pepck_c]	1	0.6	
54	pepck_c →	1	0.6	logical model output / metabolic input
55	lxra = pfk2_fbp2_g [reacs_pfk2_fbp2_g]	1	1.0	LXRα was shown to transactivate the pfkfb1 promoter via LXRE binding and to up-regulate PFK2 / FBP2 gene expression upon insulin treatment or refeeding (as demonstrated in murine liver and for HuH7 cells146)

56	nsrebp1c = pfk2_fbp2_g [reacs_pfk2_fbp2_g]	1	0.8	SREBP1 transactivates the pfkfb1 promoter via SRE binding in vivo, likely mediating the up-regulation of PFK2 / FBP2 gene expression in response to insulin or glucose stimulation (as demonstrated in livers of S. aurata80)

57	nsrebp1c = lpin1 [reacs_lpin1]	1	0.6	nSREBP1c and NF-Y cooperatively induce LPIN1 gene expression (as shown by co-transfection studies in HepG2 cells; nSREBP1 and NF-Y binding to lpin1 promoter demonstrated for LPDS-treated Huh7 cells55); NF-Y involvement omitted for model reduction purposes
58	pgc1a = lpin1 [reacs_lpin1]	1	0.8	PGC1α was shown to be required and sufficient for fasting-induced activation of LPIN1 gene expression (as shown in liver of fasted or dexamethasone-treated mice36)

Glycolysis / gluconeogenesis
59	!pka · pp2a · pfk2_fbp2_g = pfk2 [reacs_pfk2]	1	1.0	as generally accepted, PKA induces negative-regulatory PFK2-S32 phosphorylation (residue of rat liver isoform!) in response to glucagon stimulation, leading to conformational changes of the bifunctional enzyme, the latter in turn increasing FBP2 activity; PP2A (e.g. activated via glucose / X5P87) counteracts PKA action, dephosphorylating PFK2-S32 and inducing PFK2 activity (for review see [89, 90])

60	pka · !pp2a · pfk2_fbp2_g = fbp2 [reacs_fbp2]	1	1.0	
61	pfk2 →	1	1.0	logical model output / metabolic input
62	fbp2 →	1	1.0	logical model output / metabolic input
63	!pka · lpk_g = lpk [reacs_lpk]	1	0.6	PKA catalyses negative-regulatory L-PK serine phosphorylation, decreasing L-PK substrate / allosteric effector (PEP, F-1,6-P2) affinity and subsequent kinase activity
64	lpk →	1	0.6	logical model output / metabolic input
65	→ pyruvate [reacs_pyruvate]	1	1.0	logical model input / metabolic output
66	!pyruvate · pdhk2_g = pdhk2 [reacs_pdhk2]	1	1.0	pyruvate inhibits PDHK activation via direct interaction12, conceivably impeding PDHK autophosphorylation; for review see [119]

67	npkcd = pdhp2 [reacs_pdhp2]	1	0.8	insulin induces mitochondrial nPKCδ translocation, leading to nPKCδ:PDHP2 interaction and nPKCδ-catalysed PDHP2 phosphorylation / activation, in turn upregulating PDHC activity within 10 min (as shown for insulin-treated immortalised murine hepatocytes16)

68	!pdhk2 · pdhp2 = pdhc [reacs_pdhc]	1	1.0	PDHK2 catalyses the negative-regulatory phosphorylation on inactivating sites 1 and 2 (S264, S271) within the PDHC-E1α subunit, whereas dephosphorylation by PDHP2 results in PDHC reactivation (as demonstrated in vitro66; for review see [119])

69	pdhc →	1	1.0	logical model output / metabolic input
 
Glycogen synthesis / glycogenolysis
70	!pp2a · pka = pygl [reacs_pygl]	1	0.8	PHK catalyses positive-regulatory PYGL-S14 phosphorylation, resulting in “T-to-R state” transition and PYGL activation (for review see [38]); PKA-catalysed phosphorylation of the regulatory β subunits (likely inducing subsequent α subunit phosphorylation) directly increases PHK activity (as shown for rabbit skeletal muscle PHK98; for review see [13]); according to their classification, PP1 and PP2A selectively dephosphorylate PHK β and α subunits, respectively (for review see [53]); explicit PHKL involvement omitted / bypassed for model reduction purposes
71	pygl →	1	0.8	logical model output / metabolic input
72	pp2a · !gsk3 · !pka · gys2_g = gys2 [reacs_gys2]	1	0.6	PP1-GL dephosphorylates (as shown for GSK3 target sites in vitro29) and thus activates GYS (for review see [96]); generally, substrate specificities of PP1 and PP2A for GSK3 and PKA target sites have been demonstrated in skeletal muscle (for review see [53]); insulin was shown to activate glycogen synthase via GSK3 inactivation, hence preventing negative-regulatory GYS S641,  645,  549, and  653 phosphorylation (for review see [99]); additionally, insulin seems to induce phosphatases that directly dephosphorylate those residues targeted by GSK3; PKA was shown to catalyse the inhibitory phosphorylation at sites 1a/b and 2 of rabbit skeletal muscle glycogen synthase in vitro94; but: individual relevances of inhibitory kinases for hepatocellular GYS regulation have to be checked; PP1 involvement omitted for model reduction purposes
'''

special_variables = dict({'dum_cbp_or_torc2_pepck_c':'cbp OR torc2',
      'dum_cbp_or_torc2_pgc1a_g':'cbp OR torc2',
      'dum_gluc_in_or_g6p_lxra':'gluc_in OR g6p',
      'dum_lxra_or_basal_exp_insig2a':'lxra OR basal_exp',
      'dum_lxra_sp1_psrebp1c':'lxra * sp1',
      'dum_hnf4a_or_ppara_ppara_g':'hnf4a OR ppara',
      'dum_cbp_or_torc2_g6pc':'cbp OR torc2',
      'dum_cbp_or_torc2_pepck_c':'cbp OR torc2'})

lines = text.replace('∙','AND').replace('·','AND').replace('!','NOT ').splitlines()




count = 1
variables,rules = [],[]
for line in lines:
    if line=='':
        continue
    try:
        index = line.index('[')
        line = line[:index]
    except ValueError:
        pass
    linesplit = line.split('\t')
    if linesplit[0]==str(count): #next rule
        print(linesplit[1])
        count+=1
        if not '=' in linesplit[1]: #constant
            continue
        rule,var = linesplit[1].split('=')
        for special_var in special_variables.keys():
            if special_var in rule:
                rule = rule.replace(special_var,'( '+special_variables[special_var]+' )')
        try:
            index_var = variables.index(var)
            rules[index_var] += ' OR '+'( '+rule+' )'
        except:
            variables.append(var)
            rules.append('( '+rule+' )')

for first_star in ['AND','OR']:
    for second_star in ['AND','OR']:
        
        counter_star = 0
        g = open(folder+pmid+'_%s_%s.txt' % (first_star,second_star),'w')
        for var,rule in zip(variables,rules):
            if '*' in rule:
                rule = rule.replace('*',first_star if counter_star==0 else second_star)
                counter_star +=1
            g.write((var+' = '+rule+'\n').replace('  ',' '))
        g.close()        

