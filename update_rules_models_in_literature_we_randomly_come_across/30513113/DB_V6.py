'''
This version implements the final version of DB in their 2013 paper
'''


import itertools

def findNextCycle(preStart,preCig1Cdc2,preCig2Cdc2,prePuc1Cdc2,preCdc2Cdc13,preSte9,preRum1,preSlp1,preCdc2_Tyr15,preWee1Mik1,preCdc25,prePP):
    value=""

    #Start
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
'''
state='100001100100'
print "1 "+str(state)
for i in range(2,20):
    new_state=findNextCycle(int(state[0]),int(state[1]),int(state[2]),int(state[3]),int(state[4]),int(state[5]),int(state[6]),int(state[7]),int(state[8]),int(state[9]),int(state[10]),int(state[11]))
    print str(i)+" "+str(new_state)
    state=new_state
'''
print " start "+"cig1/cdc2"+"cig2/cdc2"+"puc1/cdc2"+"cdc2cdc13"+"  ste9 "+"  rum1 "+"  slp1 "+"cdc2tyr15"+"wee1/mik1"+" cdc25 "+"   pp   "
state='100001100100'
print "1"+"    "+str(state[0])+"       "+str(state[1])+"       "+str(state[2])+"       "+str(state[3])+"       "+str(state[4])+"       "+str(state[5])+"       "+str(state[6])+"       "+str(state[7])+"       "+str(state[8])+"       "+str(state[9])+"       "+str(state[10])+"       "+str(state[11])
for i in range(2,20):
    new_state=findNextCycle(int(state[0]),int(state[1]),int(state[2]),int(state[3]),int(state[4]),int(state[5]),int(state[6]),int(state[7]),int(state[8]),int(state[9]),int(state[10]),int(state[11]))
    #print str(i)+" "+str(new_state)
    print str(i)+"    "+str(new_state[0])+"       "+str(new_state[1])+"       "+str(new_state[2])+"       "+str(new_state[3])+"       "+str(new_state[4])+"       "+str(new_state[5])+"       "+str(new_state[6])+"       "+str(new_state[7])+"       "+str(new_state[8])+"       "+str(new_state[9])+"       "+str(new_state[10])+"       "+str(new_state[11])
    state=new_state
