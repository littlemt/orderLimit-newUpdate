#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 15:23:02 2022

@author: Milo Kamiya Belmont
"""

import dep.FPDMCUpdates_OrderLimit as FPC
import time
import numpy as np
import numpy.random as nrng
import matplotlib.pyplot as mpl

    

nrg=nrng.default_rng()


def full(tauMax,runTime,p,pExt,mu,k,alpha,orderMax,m=1):
    
    
    tauList=[]
    qList=np.ndarray((orderMax,5))
    
    total=sum(p)
    
    pTau=p[0]/total
    pIns=p[1]/total
    pRem=p[2]/total
    pSwap=p[3]/total
    
    countT=0
    countI=0
    countR=[]
    countS=0
    tau=tauMax
    mcTime=[]
    mcT=1
    
    orderList=[]
    runTime=runTime*3600+time.time()
    
    while time.time()<runTime:
        
        
        order=FPC.order(qList)
        orderList.append(order)
        
        
        x=nrg.uniform()
        
        
        
        if x<pTau and order==0:
            
            mcTime.append(mcT)
            tau,i=FPC.tauChange(pExt,tauMax,tau,mu,m)
            tauList.append(tau)
            countT+=i
            mcT=0
        elif pTau<=x<pTau+pIns:
            
            qList,i=FPC.insertProp(qList, tau, k, mu, alpha, order,pIns,pRem)
            countI+=i
            
            
        elif pTau+pIns<=x<pTau+pIns+pRem and order>=1:
            
            qList,i=FPC.removeProp(qList,tau, k, mu, alpha, order,pIns,pRem)
            countR.append(i)
            
            
        elif pTau+pIns+pRem<=x<1 and order>=2:
            
            qList,i=FPC.swap(qList, k, tauMax, mu, order)
            countS+=i
            
        mcT+=1
            
            
               
    mcTime.append(mcT)
    count=[countT,countI,countR,countS]
    return tauList,mcTime,qList,count,orderList
        

def zero_order(tauMax,runTime,pExt,mu,alpha,m=1):
    tauList=[]
    runTime=runTime*3600+time.time()
    count=0
    tau=0
    while time.time()<runTime:
        tau,i = FPC.changeTau(tau,tauMax,pExt,0,m)
        tauList.append(tau)
        count += i
        
    return tauList,count


def first_order(tauMax,runTime,P,pExt,mu,k,alpha,orderMax,omega=1,m=1):
    '''
    

    Parameters
    ----------
    tauMax : float
        Maximum allowed tau value allowed to be picked.
    runTime : float
        Total time the simulation should run.
    P : list or array of size 3
        List defining the probabliitys for each update to be picked each loop.
    pExt : float
        External momentum of the system.
    mu : float
        Chemicle potnetial of the system.
    k : float
        Momentun of the bare electron propogator.
    alpha : float
        coupling constant.
    orderMax : int
        maximum order allowed for the simulation.
    omega : float, optional
        frequency of the particle. The default is 1.
    m : float, optional
        mass of the particle. The default is 1.

    Returns
    -------
    tauList : array type
        List of each change in tau.
    mcTime : array type
        ammount of monte carlo time between each tau update.
    qList : array type (Max order x 5)
        the final array that contains all the phonon arcs created by the simulation.
    orderList : array type
        list of the order at each update.
    count : list 
        List of the amount of time each update was accepted in order [change tau,
        insert, -remove].

    '''
    qList=np.ndarray((orderMax,5))
    tau=FPC.changeTau(0,tauMax,qList,pExt,0,m)
    tauList=[tau[0]]
    
    total=sum(P)
    pTau=P[0]/total
    pIns=P[1]/total
    pRem=P[2]/total
    countT=0
    countI=0
    countR=0
    orderList=[]
    mcTime=[]
    mcT=1
    n=0
    #startTime=time.time()
    endTime=runTime*3600+time.time()
    
    while time.time()<endTime:
        #if time.time()-startTime
        x=nrg.uniform()
        
        #print('q',qList)
        
        if 0<=x<=pTau:
            tau,i = FPC.changeTau(tau,tauMax,qList,pExt,n,m)
            tauList.append(tau)
            countT += i
            mcTime.append(mcT)
            mcT=0
        elif pTau<x<=pTau+pIns and n<orderMax:
            qList,i=FPC.insertArc(qList,pExt,tauMax,orderMax,omega,m,n)
            countI+=i
            n+=i
        elif pTau+pIns<x<=1 and n>=1:
            qList,i=FPC.removeArc(qList,omega,tauMax,orderMax,m,pExt,n)
            countR+=i
            n+=i
        orderList.append(n)  
        mcT+=1
        
    mcTime.append(mcT)
    count=[countT,countI,countR]
    return tauList,mcTime,qList,orderList,count

#check out end of second talk
#check how often a update is being rejected 
#check to see if first order is working 


def data_Unravel(qList):
    length=len(qList)
    
    table=np.ndarray((length,3))
    for i in range(length):
        q,tau1,tau2=qList[i]
        table[i]=[q,tau1,tau2]
        
    return table

def countZero(orderList):
    return np.count_nonzero(orderList==0)

def calc(tauList,mctList,noBin,tauMax,k,mu,zeroOrder,m=1):
    
    
    len1=len(tauList)
    len2=len(mctList)
    
    if len1!=len2:
        print('List length not equal')
        return
    
    timeArray=np.ndarray((len1,2))
    timeArray[:,0]=tauList
    timeArray[:,1]=mctList
    
    histBin=np.ndarray((noBin,3))
    histBin[:,0]=np.linspace(0, tauMax,noBin)
    
    deltaTau=tauMax/noBin
    
    
    for i in range(noBin):
        count=0
        for j in timeArray:
            if i*deltaTau<=j[0]<(i+1)*deltaTau:
                count+=j[1]
        
        histBin[i,1]=count
            
     
    
    epsK=k**2/(2*m)
    
    integral=1/(epsK-mu)*(np.exp(-(epsK-mu)*tauMax)-1)
    
    for i in range(noBin):
        histBin[i,2]=-histBin[i,1]*integral/deltaTau/-zeroOrder

    return histBin        
        
        
    
        

        
#add extend 



#if i limit the order then i can use numpy array
#this would require the order function to be changed would need to store 
#the order as a seperate variable or can remake the list every time change. 
#can make paralell by using multipul seeds then running the simulation in parallel it looks like

        
        