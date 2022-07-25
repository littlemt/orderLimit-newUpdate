#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 19:56:57 2022

@author: Milo Koyama Belmont
"""

import numpy as np

nrand=np.random.default_rng()

'''
This code is going to be built from the ground up with an order limit in mind.
It will also have the updates from https://youtu.be/-Zr4xYJWHbs

insert, remove
change time
swap


'''
    

def changeTau(t,tauMax,pExt,order,m):
    if order!=0:
        eps=pExt**2/(2*m)
        
        R=nrand.uniform()
        tauNew=t-np.log(R)/eps
        
        return tauNew,1
    else:
        eps=pExt**2/(2*m)
        
        R=nrand.uniform()
        tauNew=-np.log(R)/eps
        return tauNew,1
    
    

def insertArc(qList,p,tMax,orderMax,omega,m,n):
    
    #n=order(qList)
    
    tList=qList[:n,0:2]
    tList=tList.flatten()
    #tList=np.trim_zeros(tList)
    tList=np.append(tList, [0,tMax])
    
    tList=np.sort(tList,axis=None)
    
    index=nrand.integers(0,tList.size-1)
    
    
    tauOne=tList[index]
    tauOneP=tList[index+1]
    tauTwo=nrand.uniform(tauOne,tauOneP)
    tauTwoP=tauTwo+nrand.exponential(1/omega,None)
    if tauTwoP>tMax:
        return qList,0
    
    
    dTa=tauTwoP-tauTwo
    
    sigma=(m/(tauTwoP-tauTwo))**.5
    qTwo=nrand.normal(0,sigma)
    
    r=R_insert(omega, m, p, tauTwo, tauTwoP, qTwo, n, dTa)
    
    x=nrand.uniform()
    
    if x<r:
        qList[n]=np.array([tauTwo,tauTwoP,qTwo,0,0])
         
        #need to check if i am splicing right
        return qList,1
    else:
        return qList,0
        
    
    
    
def removeArc(qList,omega,tMax,orderMax,m,p,n):
    #n=order(qList)
    
    i=nrand.integers(0,n)
    [tauTwo,tauTwoP,q1,q2,q3]=qList[i]
    qTwo=np.linalg.norm([q1,q2,q3])
    
    #issue is that i am getting 
    tList=qList[:n,0:2]
    tList=tList.flatten()
    tList=np.append(tList, [0,tMax])
    tList=np.sort(tList,axis=None)
    
    #might have to make it order -1
    #how do i calc dTa, going to pick arc that covers tauTwo 
    j=np.where(tList==tauTwo)[0]
    #print(tList)
    
    #somehow my tauTwo is picking zero should not be an issue with insert 
    #therefor must be a problem with remove cant figure out where 
    
    if tList[j+1]==tauTwoP:
        tauOne=tList[j-1]
        tauOneP=tList[j+2]
    else:
        tauOne=tList[j-1]
        tauOneP=tList[j+1]
        
    dTa=tauOneP-tauOne
    r=R_insert(omega, m, p, tauTwo, tauTwoP, qTwo, n, dTa)**-1
    
    x=nrand.uniform()
    #print(i,n)
    dum1=qList[i+1:n+1]
    dum2=qList[i:n]
    #print(np.shape(dum1),np.shape(dum2))
    if x<r:
        print('a',i)
        if i==n-1:
            #this part says if the arc to remove is the same as the max order 
            #then just remove the last arc from the list 
            qList[n-1]=np.zeros(5)
        else:
            #takes the arcs on the list from after the arc picked and places them back one spot
            
            dummy=qList[i+1:orderMax]
            qList[i:orderMax-1]=dummy
            qList[orderMax-1]=np.zeros(5)
        
        return qList,-1
    else:
        print('r')
        return qList,0
        
        
        
def swap (qList,order):
    
    a=nrand.integers(0,order)
    [b,c]=nrand.integers(0,2,size=2)
    
    arcOne=qList[a]
    
    tauList=qList[:order,:2]
    tauList=tauList.flatten()
    tauList=np.sort(tauList)
    
    
    if arcOne(b)==tauList[0]:
        i=1
    elif arcOne[b]==tauList[-1]:
        i=tauList.size()-1
    else:
        for i in range (1,tauList.size()-1):
            if arcOne[b]==tauList[i]:
                if c==0:
                    i=i-1
                else:
                    i=i+1
            break
    
    x=nrand.uniform()
    
    ind=np.where(qList[:order,:2]==tauList[i])
    arcTwo=qList[ind[0]]
    
    r=np.exp(-abs(arcOne[b]-arcTwo[ind[1]]))
            #eps(p')-esp(p)\pm omega_a\pm omega_b
            #what does the plusminus do? 
    if x<r:
        qList[a,b]=arcTwo[ind[1]]
        qList[a,2:5]=arcTwo[2:5]
        qList[ind]=arcOne[b]
        qList[ind[0],2:5]=arcOne[2:5]
        return qList,1
    else:
        return qList,0
        
    '''
    questions
    does this need to be done with only the initial point or is it anypoint on an interval
    and i can swap the endpoints even if they are on opposite ends. 
    
    
    
    '''
    
    
    
    
    
    
def R_insert(omega,m,p,tauTwo,tauTwoP,qTwo,n,dTa):
    #D_n+1 /D_n /Omega *(2n+1)deltaTau/
    #Omega=
    #D_n+1/D_n=
    
    D=np.exp(-omega*(tauTwoP-tauTwo))*np.exp(-(qTwo**2-2*p*qTwo)/(2*m)*(tauTwoP-tauTwo))
    qNot=np.sqrt(2*m*omega)
    Omega=1/(4*np.pi)/qNot*omega*np.exp(-omega*(1+qTwo/qNot)**2*(tauTwoP-tauTwo))
    return D/Omega*(2*n+1)*(2*n+1)*dTa/(n+1)

    
'''    
def order(qList):
    x=np.count_nonzero(qList[0])
    return x
    '''
    
            
    
#

        