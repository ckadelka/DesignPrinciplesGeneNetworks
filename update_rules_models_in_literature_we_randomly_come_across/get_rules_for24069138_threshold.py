#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:06:02 2019

@author: ckadelka
"""




var = ['Start','Cig1_Cdc2','Cig2_Cdc2','Puc1_Cdc2','Rum1','Cdc2_Cdc13','Ste9','PP','Cdc25','Slp1','Wee1','Cdc2_Tyr15']
thresholds = [0,0,0,0,0,-.5,0,0,0,1,0,0.5]

n = len(var)

A = np.zeros((n,n))
A[0,0] = 1
A[1,1] = 1
A[2,2] = 1
A[3,3] = 1
A[7,7] = 1
A[9,9] =1

A[0,1] =1
A[0,2] =1
A[0,3] =1

A[7,4] =1
A[9,7] =1
A[9,5] =1
A[11,9] =1
A[7,10] =1
A[7,6] =1
A[5,8] =1
A[8,11] =1

A[1,4] = -1
A[1,6] = -1
A[2,4] = -1
A[2,6] = -1
A[3,4] = -1
A[3,6] = -1
A[4,5] = -1
A[5,4] = -1
A[5,6] = -1
A[6,5] = -1
A[9,5] = -1
A[5,10] = -1
A[10,11] = -1

def combine(m, n):
    a = len(m)
    c = ''
    count = 0
    for i in range(a): 
        if(m[i] == n[i]):
            c += m[i]
        elif(m[i] != n[i]):
            c += '-'
            count += 1

    if(count > 1): 
        return None
    else:            
        return c


def find_prime_implicants(data):
    newList = list(data)
    size = len(newList)
    IM = []
    im = []
    im2 = []
    mark = [0]*size
    m = 0
    for i in range(size):
        for j in range(i+1, size):
            c = combine( str(newList[i]), str(newList[j]) )
            if c != None:
                im.append(str(c))
                mark[i] = 1
                mark[j] = 1
            else:
                continue

    mark2 = [0]*len(im)
    for p in range(len(im)):
        for n in range(p+1, len(im)):
            if( p != n and mark2[n] == 0):
                if( im[p] == im[n]):
                    mark2[n] = 1


    for r in range(len(im)):
        if(mark2[r] == 0):
            im2.append(im[r])

    for q in range(size):
        if( mark[q] == 0 ):
            IM.append( str(newList[q]) )
            m = m+1

    if(m == size or size == 1):
        return IM
    else:
        return IM + find_prime_implicants(im2)



F = []
for i in range(n):
    ids_regulators = A[:,i]!=0
    ids_regulators[i] = True#threshold function so the variable itself always matters if there's a tie
    n_regulators = sum(ids_regulators)
    names_regulators = np.array(var)[np.array(fbn.find_all_indices(ids_regulators,True))]
    id_i = sum(ids_regulators[:i])
    bool_list = list(map(np.array,list(itertools.product([0, 1], repeat = n_regulators))))
    what_leads_to_1 = []
    
    for el in bool_list:
        AX = np.dot(A[ids_regulators,i],el)
        if AX>thresholds[i] or AX==thresholds[i] and el[id_i]==1:
            what_leads_to_1.append( ''.join(np.array(el,dtype='str_')) )
    
    dummy = find_prime_implicants(what_leads_to_1)
    AND_terms = []
    for minterm in dummy:
        dummy2 = []
        for j in range(n_regulators):
            if minterm[j]=='1':
                dummy2.append(names_regulators[j])
            elif minterm[j]=='0':
                dummy2.append('NOT '+names_regulators[j])
        AND_terms.append(' AND '.join(dummy2))
        if len(dummy2)>1:
            AND_terms[-1] = '('+AND_terms[-1]+')'
    F.append(' OR '.join(AND_terms))


g = open('update_rules_models_in_literature_we_randomly_come_across/24069138_threshold.txt','w')
for i in range(n):
    g.write(var[i]+' = '+F[i]+'\n')
g.close()

