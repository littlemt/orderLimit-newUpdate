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


def first_order(tauMax,runTime,P,pExt,mu,alpha,orderMax,omega=1,m=1):
    '''
    

    Parameters
    ----------
    tauMax : float
        Maximum allowed tau value allowed to be picked.
    runTime : float
        Total time the simulation should run.
    P : list or array of size 4
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
    qList=np.zeros((orderMax,3))
    mList=np.zeros((orderMax*2+2,2))
    
    tau=FPC.changeTau(0,tauMax,mList,pExt,0,mu,m)
    tauList=[tau[0]]
    mList[0:2,0]=[0,tau[0]]
    mList[0,0]=pExt
    
    total=sum(P)
    pTau=P[0]/total
    pIns=sum(P[:2])/total
    pRem=sum(P[:3])/total
    pSwap=sum(P[:4])/total
    pExt=sum(P[:5])/total
    print(pTau,pIns,pRem,pSwap)
    countT=0
    countI=0
    countR=0
    countS=0
    countE=0
    orderList=[]
    mcTime=[]
    mcT=1
    n=0
    #startTime=time.time()
    endTime=runTime*3600+time.time()
    
    while time.time()<endTime:
        #if time.time()-startTime
        x=nrg.uniform()
        #print(pTau<x<pIns,x)
        #print('q',qList)
        #print('m',mList)
        #print(n)
        if 0<=x<pTau and n==0:
            #print('tau')
            tau,i = FPC.changeTau(tau,tauMax,mList,pExt,n,mu,m)
            tauList.append(tau)
            countT += i
            mList[2*n+1,0]=tau
            #print('tau')
            mcTime.append(mcT)
            mcT=0
        elif pTau<x<=pIns and n<orderMax:
            #print('ins')
            qList,mList,i=FPC.insertArc(qList,mList,tau,orderMax,omega,m,n,pIns,pRem,alpha,mu)
            countI+=i
            n+=i
            
        elif pIns<x<=pRem and n>=1:
            #print('rem')
            qList,mList,i=FPC.removeArc(qList,mList,orderMax,omega,m,n,mu,pRem,pIns,alpha)
            countR+=i
            n+=i
        elif pRem<=x<pSwap and n>=2:
            
            qList,mList,i=FPC.swap(qList, n)
            countS+=i
            
        elif pSwap<=x<pExt and n<=1:
            tau,i = FPC.changeTau(tau,tauMax,mList,pExt,n,mu,m)
            tauList.append(tau)
            countT += i
            mList[2*n+1,0]=tau
            #print('tau')
            mcTime.append(mcT)
            mcT=0
        
        orderList.append(n)  
        mcT+=1
        
    mcTime.append(mcT)
    count=[countT,countI,countR]
    return tauList,mcTime,qList,orderList,count,mList

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

def calc(tauList,mctList,noBin,tauMax,k,mu,zeroOrder,thermal,skip,m=1):
    '''
    

    Parameters
    ----------
    tauList : TYPE
        DESCRIPTION.
    mctList : TYPE
        DESCRIPTION.
    noBin : int
        number of bins.
    tauMax : float
        maximum allowed external time.
    k : float
        DESCRIPTION.
    mu : float
        DESCRIPTION.
    zeroOrder : int
        number of times in zero order.
    m : float, optional
        mass. The default is 1.

    Returns
    -------
    histBin : array type (noBinx3)
        returns an array of bins with collum one being the bins, collum two is
        the count of the amount of updates that happen in the external time in that bin,
        and collum three is the calculated value for that bin

    '''
    
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
        for j in range(thermal,len(timeArray),skip):
            if i*deltaTau<=j[0]<(i+1)*deltaTau:
                count+=j[1]
        
        histBin[i,1]=count
            
     
    
    epsK=k**2/(2*m)
    
    integral=1/(epsK-mu)*(np.exp(-(epsK-mu)*tauMax)-1)
    
    for i in range(noBin):
        #just realized this is implimented wrong i think
        histBin[i,2]=-histBin[i,1]*integral/deltaTau/-zeroOrder

    return histBin        
        
        
    
        

        
#add extend 



#if i limit the order then i can use numpy array
#this would require the order function to be changed would need to store 
#the order as a seperate variable or can remake the list every time change. 
#can make paralell by using multipul seeds then running the simulation in parallel it looks like

        
        