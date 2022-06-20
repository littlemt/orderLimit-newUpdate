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
        
    
    
    
def removeArc(qList,omega,tMax,m,p,n):
    #n=order(qList)
    
    i=nrand.integers(0,n)
    [tauTwo,tauTwoP,q1,q2,q3]=qList[i]
    qTwo=np.linalg.norm([q1,q2,q3])
    
    
    tList=qList[:n,0:2]
    tList=tList.flatten()
    #tList=np.trim_zeros(tList)
    tList=np.append(tList, [0,tMax])
    tList=np.sort(tList,axis=None)
    
    #might have to make it order -1
    #how do i calc dTa, going to pick arc that covers tauTwo 
    j=np.where(tList==tauTwo)[0]
    print(tList)
    
    #i know what is happening here dont know how to fix it. maybe slicing wrong?
    print(tList[j+1],j+1)
    if tList[j+1]==tauTwoP:
        tauOne=tList[j-1]
        tauOneP=tList[j+2]
    else:
        tauOne=tList[j-1]
        tauOneP=tList[j+1]
        
    dTa=tauOneP-tauOne
    r=R_insert(omega, m, p, tauTwo, tauTwoP, qTwo, n, dTa)**-1
    
    x=nrand.uniform()
    print(i,n)
    dum1=qList[i+1:n+1]
    dum2=qList[i:n]
    print(np.shape(dum1),np.shape(dum2))
    if x<1:
        print('a')
        if i==n-1:
            qList[n-1]=np.zeros(5)
        else:
            dummy=qList[i+1:n+1]
            qList[n-1]=np.zeros(5)
            qList[i:n]=dummy
        
        return qList,-1
    else:
        print('r')
        return qList,0
        
        
        
    
    
    
    
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
    
            
    
    #need to create list of nodes
    
    #new change tau same as old but can do at not 0 order 1:03:40
    #new insert remove 
    #max order 2 
    #add extend

        