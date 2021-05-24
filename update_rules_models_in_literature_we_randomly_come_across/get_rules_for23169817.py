#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""

HIGH_DAMAGE=False

if HIGH_DAMAGE==False:
    #low DNA damage rules
    a = '''ATM*																										
    e2f1	and	not	cycg	and	not	wip1	and	atm(t)																		
    p53*																										
    not	atm	and	not	mdm2	and	not	mdmx																			
    atm	and	not	mdm2																							
    atm	and	mdm2	and	not	mdmx	and	p53(t)																			
    wip1*																										
    p53																										
    cycg*																										
    p53																										
    pten*																										
    p53																										
    p21*																										
    p53	and	not	akt																							
    p53	and	akt	and	not	mdm2																					
    p53	and	akt	and	mdm2	and	p21(t)																				
    AKT*																										
    not	pten																									
    cycE*																										
    not	p21																									
    rb*																										
    not	atm	and	not	cyce	and	not	mdm2	and	not	casp	and	rb(t)														
    atm	and	not	cyce	and	not	casp																				
    atm	and	cyce	and	not	mdm2	and	not	casp																		
    atm	and	cyce	and	mdm2	and	not	casp	and	rb(t)																	
    e2f1*																										
    not	arf	and	not	rb																						
    not	arf	and	rb	and	not	atm	and	not	mdm2	and	e2f1(t)															
    not	arf	and	rb	and	not	atm	and	mdm2																		
    not	arf	and	rb	and	atm																					
    arf	and	not	rb	and	atm	and	mdm2	and	e2f1(t)																	
    arf*																										
    not	p53	and	not	wip1	and	e2f1	and	arf(t)																		
    bcl2*																										
    not	akt	and	not	p53	and	not	casp	and	bcl2(t)																	
    akt	and	not	p53																							
    akt	and	p53	and	not	casp	and	bcl2(t)																			
    bax*																										
    p53	and	not	bcl2																							
    p53	and	bcl2	and	bax(t)																						
    casp*																										
    not	akt	and	not	bcl2	and	not	p21	and	not	bax	and	e2f1	and	casp(t)												
    not	akt	and	not	p21	and	bax																				
    not	akt	and	not	bcl2	and	p21	and	bax	and	not	e2f1	and	not	casp(t)												
    not	akt	and	p21	and	bax	and	e2f1																			
    akt	and	not	bcl2	and	not	p21	and	bax																		
    akt	and	not	bcl2	and	p21	and	bax	and	e2f1																	
    akt	and	bcl2	and	not	p21	and	bax	and	not	e2f1	and	casp(t)														
    akt	and	bcl2	and	not	p21	and	bax	and	e2f1																	
    akt	and	bcl2	and	p21	and	bax	and	e2f1	and	casp(t)																
    mdmx*																										
    not	atm	and	not	mdm2	and	not	arf	and	not	akt																
    not	mdm2	and	not	arf	and	akt																				
    not	atm	and	not	mdm2	and	arf	and	akt	and	wip1																
    not	atm	and	mdm2	and	not	arf	and	akt	and	wip1	and	mdmx(t)														
    atm	and	not	mdm2	and	not	arf	and	not	akt	and	wip1	and	mdmx(t)													
    mdm2*																										
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1)				
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)	
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	mdmx												
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53							
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1				
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53							
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1)				
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx										
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg													
    not	arf	and	not	atm	and	not	cyce	and	not	rb																
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1)						
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53									
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1						
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53									
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1)						
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1							
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx												
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1				
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)	
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	mdmx												
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53							
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1)				
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53							
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1				
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx										
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg													
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1						
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53									
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1)						
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53									
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1)						
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1)							
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx												
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1)				
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)	
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx									
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53							
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1)				
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx										
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53							
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1)				
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)	
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx										
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg													
    not	arf	and	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)						
    not	arf	and	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)						
    not	arf	and	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53	and	wip1)									
    not	arf	and	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)						
    not	arf	and	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)						
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1						
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53									
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1)						
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53									
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1						
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1)							
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx												
    arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)				
    arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	wip1)							
    arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	wip1)							
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)					
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53										
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	wip1)							
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)				
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)		
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53							
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1)				
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)	
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	wip1)					
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)		
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)		
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx										
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1)					
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)		
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)		
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx										
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53								
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1					
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)		
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	mdmx											
    arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)							
    arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)				
    arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	wip1)							
    arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	wip1)							
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)					
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	p53										
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	wip1							
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)				
    arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)							
    arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)				
    arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	wip1)							
    arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	wip1)							
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)				
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)				
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)					
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	p53										
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	wip1							
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)				
    arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)							'''
else:
    #high DNA damage rules
    a = '''ATM*																									
    not	e2f1	and	not	cycg	and	not	wip1																		
    not	e2f1	and	not	cycg	and	wip1	and	atm(t)																	
    not	e2f1	and	cycg	and	not	wip1	and	atm(t)																	
    e2f1	and	not	cycg	and	not	wip1																			
    e2f1	and	not	cycg	and	wip1																				
    e2f1	and	cycg	and	not	wip1																				
    p53*																									
    not	atm	and	not	mdm2	and	not	mdmx																		
    atm	and	not	mdm2																						
    atm	and	mdm2	and	not	mdmx	and	p53(t)																		
    wip1*																									
    p53																									
    cycg*																									
    p53																									
    pten*																									
    p53																									
    p21*																									
    p53	and	not	akt																						
    p53	and	akt	and	not	mdm2																				
    p53	and	akt	and	mdm2	and	p21(t)																			
    AKT*																									
    not	pten																								
    cycE*																									
    not	p21																								
    rb*																									
    not	atm	and	not	cyce	and	not	mdm2	and	not	casp	and	rb(t)													
    atm	and	not	cyce	and	not	casp																			
    atm	and	cyce	and	not	mdm2	and	not	casp																	
    atm	and	cyce	and	mdm2	and	not	casp	and	rb(t)																
    e2f1*																									
    not	arf	and	not	rb																					
    not	arf	and	rb	and	not	atm	and	not	mdm2	and	e2f1(t)														
    not	arf	and	rb	and	not	atm	and	mdm2																	
    not	arf	and	rb	and	atm																				
    arf	and	not	rb	and	atm	and	mdm2	and	e2f1(t)																
    arf*																									
    not	p53	and	not	wip1	and	not	e2f1	and	arf(t)																
    not	p53	and	not	wip1	and	e2f1																			
    not	p53	and	wip1	and	e2f1	and	arf(t)																		
    bcl2*																									
    not	akt	and	not	p53	and	not	casp	and	bcl2(t)																
    akt	and	not	p53																						
    akt	and	p53	and	not	casp	and	bcl2(t)																		
    bax*																									
    p53	and	not	bcl2																						
    p53	and	bcl2	and	bax(t)																					
    casp*																									
    not	akt	and	not	bcl2	and	not	p21	and	not	bax	and	e2f1	and	casp(t)											
    not	akt	and	not	p21	and	bax																			
    not	akt	and	not	bcl2	and	p21	and	bax	and	not	e2f1	and	not	casp(t)											
    not	akt	and	p21	and	bax	and	e2f1																		
    akt	and	not	bcl2	and	not	p21	and	bax																	
    akt	and	not	bcl2	and	p21	and	bax	and	e2f1																
    akt	and	bcl2	and	not	p21	and	bax	and	not	e2f1	and	casp(t)													
    akt	and	bcl2	and	not	p21	and	bax	and	e2f1																
    akt	and	bcl2	and	p21	and	bax	and	e2f1	and	casp(t)															
    mdmx*																									
    not	atm	and	not	mdm2	and	not	arf	and	not	akt															
    not	mdm2	and	not	arf	and	akt																			
    not	atm	and	not	mdm2	and	arf	and	akt	and	wip1															
    not	atm	and	mdm2	and	not	arf	and	akt	and	wip1	and	mdmx(t)													
    atm	and	not	mdm2	and	not	arf	and	not	akt	and	wip1	and	mdmx(t)												
    mdm2*																									
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1)			
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	mdmx											
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53						
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1			
    not	arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53						
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1)			
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx									
    not	arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg												
    not	arf	and	not	atm	and	not	cyce	and	not	rb															
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1)					
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)		
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)			
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53								
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1					
    not	arf	and	not	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)			
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53								
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1)					
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1						
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    not	arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx											
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1			
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	mdmx											
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53						
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1)			
    not	arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53						
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1			
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx									
    not	arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg												
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1					
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)		
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)			
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53								
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1)					
    not	arf	and	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)			
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53								
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1)					
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1)						
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    not	arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx											
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1)			
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx								
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53						
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1)			
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx									
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53						
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1)			
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx									
    not	arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg												
    not	arf	and	atm	and	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)					
    not	arf	and	atm	and	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)					
    not	arf	and	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53	and	wip1)								
    not	arf	and	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)					
    not	arf	and	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)					
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1					
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)		
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)			
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53								
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1)					
    not	arf	and	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)			
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53								
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1					
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)		
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1)						
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    not	arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx											
    arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)			
    arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	wip1)						
    arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    arf	and	not	atm	and	not	cyce	and	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	wip1)						
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53									
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	wip1)						
    arf	and	not	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)	
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53						
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1)			
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	wip1)				
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)	
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)	
    arf	and	not	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx									
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	wip1)				
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	p53	and	not	wip1	and	mdm2(t)	
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	not	mdmx	and	not	p53	and	wip1	and	mdm2(t)	
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx									
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53							
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	wip1				
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)	
    arf	and	not	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	mdmx										
    arf	and	not	atm	and	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)						
    arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)			
    arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	wip1)						
    arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    arf	and	not	atm	and	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	wip1)						
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	p53									
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	wip1						
    arf	and	not	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    arf	and	atm	and	not	cyce	and	rb	and	akt	and	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)						
    arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	not	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)			
    arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	wip1)						
    arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    arf	and	atm	and	not	cyce	and	not	rb	and	not	akt	and	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	wip1)						
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	p53	and	not	wip1	and	mdm2(t)			
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	not	cycg	and	mdmx	and	not	p53	and	wip1	and	mdm2(t)			
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	not	mdmx	and	p53	and	wip1	and	mdm2(t)				
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	p53									
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	wip1						
    arf	and	atm	and	not	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	not	p53	and	not	wip1	and	mdm2(t)			
    arf	and	atm	and	cyce	and	not	rb	and	akt	and	cycg	and	mdmx	and	p53	and	wip1	and	mdm2(t)'''
    
#results in an error because some wip still have a ) at the end --> remove

aa = []
for line in a.split('\n'):
    if '*' in line:
        aa.append(line.split('*')[0].lower() + ' = ')
    else:
        if aa[-1][-2] != '=':
            aa[-1] += ' OR '
        aa[-1] += '('+line.replace('\t',' ').replace(' and ',' AND ').replace(' not ',' NOT ').replace(' or ',' OR ').replace('(t)','').replace(')','')+')'



g = open('update_rules_models_in_literature_we_randomly_come_across/23169817_%s_dna_damage.txt' % ('high' if HIGH_DAMAGE else 'low'),'w')
for line in aa:
    g.write(line.replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ')+'\n')
g.close()

