'''
This version implements the final version of DB in their 2013 paper along with our nodes
AND editing the cdc10 role from direct activation of Cig2 to indirect activation
by activating START and thus starting a new cycle, and defining its role in activating Wee1Mik1 only to G1/S but specific conditions
and roles for our network effects based on the phase
with dephosphorylation effect at mitotic exit
'''


import itertools

def findNextCycle(preStart,preCig1Cdc2,preCig2Cdc2,prePuc1Cdc2,preCdc2Cdc13,preSte9,preRum1,preSlp1,preCdc2_Tyr15,preWee1Mik1,preCdc25,prePP, preSep1, preFkh2, preAtf1, preCdc10):
    value=""
    Next_is_G1S=0
    if preStart==0 and preCig1Cdc2==1 and preCig2Cdc2==1 and prePuc1Cdc2==1 and preCdc2Cdc13==0 and preSte9==1 and preRum1==1 and preSlp1==0 and preCdc2_Tyr15==0 and  preWee1Mik1==1 and preCdc25==0 and prePP==0  :
        Next_is_G1S=1
    Next_is_G2M=0
    if preStart==0 and preCig1Cdc2==0 and preCig2Cdc2==0 and prePuc1Cdc2==0 and preCdc2Cdc13==1 and preSte9==0 and preRum1==0 and preSlp1==0 and preCdc2_Tyr15==0 and  preWee1Mik1==0 and preCdc25==1 and prePP==0  :
        Next_is_G2M=1
    Next_is_G1=0
    if preStart==0 and preCig1Cdc2==0 and preCig2Cdc2==0 and prePuc1Cdc2==0 and preCdc2Cdc13==0 and preSte9==1 and preRum1==1 and preSlp1==0 and preCdc2_Tyr15==1 and  preWee1Mik1==1 and preCdc25==0 and prePP==1  :
        Next_is_G1=1
    Next_is_M=0
    if preStart==0 and preCig1Cdc2==0 and preCig2Cdc2==0 and prePuc1Cdc2==0 and preCdc2Cdc13==1 and preSte9==0 and preRum1==0 and preSlp1==1 and preCdc2_Tyr15==1 and  preWee1Mik1==0 and preCdc25==1 and prePP==0  :
        Next_is_M=1
    Next_is_G2=0 # used instead of S phase to indicate the Fkh2 dephopshorylation
    if preWee1Mik1==1 and preCdc2Cdc13==1:
        Next_is_G2=1
        
    #print G1, G1S, G2M
    #Start
    if preCdc2_Tyr15==1 and preWee1Mik1==1: # to mark the cell cycle is over so that the cdc10 affects start
        res=preCdc10
    else:
        res=-1*preStart
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preStart)
    #Cig1/Cdc2
    res=1*preStart-1*preCig1Cdc2
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preCig1Cdc2)
    #Cig2/Cdc2
    res=1*preStart-1*preCig2Cdc2
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preCig2Cdc2)
    #Puc1/Cdc2
    res=1*preStart-1*prePuc1Cdc2
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(prePuc1Cdc2)
    #Cdc2Cdc13
    if Next_is_G2M==1:
        res=-1*preRum1-1*preSte9-1*preSlp1+1*preAtf1
    else:
        res=-1*preRum1-1*preSte9-1*preSlp1
    if res<-0.5:
        value+='0'
    elif res>-0.5:
        value+='1'
    else:
        value+=str(preCdc2Cdc13)
    #Ste9
    res=-1*preCig1Cdc2-1*preCig2Cdc2-1*prePuc1Cdc2-1*preCdc2Cdc13+1*prePP
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preSte9)
    #Rum1
    if Next_is_M==1:
        res=-1*preCig1Cdc2-1*preCig2Cdc2-1*prePuc1Cdc2-1*preCdc2Cdc13+1*prePP+1*preSep1+1*preFkh2
    else:
        res=-1*preCig1Cdc2-1*preCig2Cdc2-1*prePuc1Cdc2-1*preCdc2Cdc13+1*prePP
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preRum1)
    #Slp1
    res=1*preCdc2Cdc13+1*preCdc2_Tyr15-1*preSlp1
    if res<1:
        value+='0'
    elif res>1:
        value+='1'
    else:
        value+=str(preSlp1)
    #Cdc2_Tyr15
    res=-1*preWee1Mik1+1*preCdc25
    if res<0.5:
        value+='0'
    elif res>0.5:
        value+='1'
    else:
        value+=str(preCdc2_Tyr15)
    #Wee1Mik1
    if Next_is_G1S==1:
        res=-1*preCdc2Cdc13+1*prePP+1*preCdc10
    else:
        res=-1*preCdc2Cdc13+1*prePP
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preWee1Mik1)
    #cdc25
    res=1*preCdc2Cdc13-1*prePP
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preCdc25)
    #PP
    res=1*preSlp1-1*prePP
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(prePP)
    #Sep1
    if Next_is_G2M==1: # it is phosphorylated by cdc2/cdc13 in G2/M
        res=1*preCdc2Cdc13
    elif Next_is_G1==1:# turn off at mitotic exit
        res=-1
    else:
        res=0
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preSep1)
    #Fkh2
    if Next_is_G1S==1: # it is phosphorylated by cdc2/cdc13 in G1/S
        res=-1
    elif Next_is_G2==1: # turn off at S phase (here used at G2 phase as the S phase is not separated from G2)
        res=1
    else:
        res=0
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preFkh2)
    #Atf1
    if Next_is_G2M==1: # it is phosphorylated by cdc2/cdc13 in G2/M
        res=1*preCdc2Cdc13
    elif Next_is_G1==1: # turn off at mitotic exit
        res=-1
    else:
        res=0
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preAtf1)
    #Cdc10
    res=1 # cdc10 is always phosphorylated and no changes on its level
    if res<0:
        value+='0'
    elif res>0:
        value+='1'
    else:
        value+=str(preCdc10)
        
    return value

def findAttractors():
    poss=[]
    for i in itertools.product(['0','1'],repeat=7):
        j=''.join(i)
        states=[]
        pre=j
        while pre not in states:
            states.append(pre)
            new=findNextCycle(pre[6],pre[5],pre[4],pre[3],pre[2],pre[1],pre[0])
            pre=new
        states.append(pre)
        decstates=[]
        for k in states:
            decstates.append(int(k,2))
        poss.append(decstates)
    return poss
    

#poss=findAttractors()
#for i in poss:
#    print(i)
print " start "+"cig1/cdc2"+"cig2/cdc2"+"puc1/cdc2"+"cdc2cdc13"+"  ste9 "+"  rum1 "+"  slp1 "+"cdc2tyr15"+"wee1/mik1"+" cdc25 "+"   pp   "+"  sep1 "+"  fkh2 "+"   atf1"+"  cdc10"
state='1000011001000101'
print "1"+"    "+str(state[0])+"       "+str(state[1])+"       "+str(state[2])+"       "+str(state[3])+"       "+str(state[4])+"       "+str(state[5])+"       "+str(state[6])+"       "+str(state[7])+"       "+str(state[8])+"       "+str(state[9])+"       "+str(state[10])+"       "+str(state[11])+"       "+str(state[12])+"       "+str(state[13])+"      "+str(state[14])+"     "+str(state[15])
for i in range(2,22):
    new_state=findNextCycle(int(state[0]),int(state[1]),int(state[2]),int(state[3]),int(state[4]),int(state[5]),int(state[6]),int(state[7]),int(state[8]),int(state[9]),int(state[10]),int(state[11]),int(state[12]),int(state[13]),int(state[14]),int(state[15]))
    #print str(i)+" "+str(new_state)
    print str(i)+"    "+str(new_state[0])+"       "+str(new_state[1])+"       "+str(new_state[2])+"       "+str(new_state[3])+"       "+str(new_state[4])+"       "+str(new_state[5])+"       "+str(new_state[6])+"       "+str(new_state[7])+"       "+str(new_state[8])+"       "+str(new_state[9])+"       "+str(new_state[10])+"       "+str(new_state[11])+"       "+str(new_state[12])+"       "+str(new_state[13])+"      "+str(new_state[14])+"     "+str(new_state[15])
    state=new_state
