#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 15:23:02 2022

@author: Milo Kamiya Belmont
"""

import dep.FPDMCUpdates_OrderLimit as FPC
import numpy as np
import numpy.random as nrng

from scipy.special import erf

import configparser
config=configparser.ConfigParser()

    




        




def first_order(tauMax,runTime,P,pExt,mu,alpha,orderMax,thermal,step,seed,mcTMax=-1,bins=100,omega=1,m=1,debug=0):
    '''
    data=first_order(20,1,[100,10,10,.1,1,0],0,-6,5,2,1000000,500,42)

    Parameters
    ----------
    tauMax : float
        Maximum allowed tau value allowed to be picked.
    runTime : float
        Total time the simulation should run.
    P : list or array of size 6
        List defining the probabliitys for each update to be picked each loop.
    pExt : float
        External momentum of the system.
    mu : float
        Chemicle potnetial of the system.
    alpha : float
        coupling constant.
    orderMax : int
        maximum order allowed for the simulation.
    mcTMax : int
        resets the loop with after a specified value 
    thermal : int
        
    omega : float, optional
        frequency of the particle. The default is 1.
    m : float, optional
        mass of the particle. The default is 1.
    debug : boolian
        if 1 then the loop is in debug mode and reports a ton of data if 0 then
        the loop only reports the nessisary data
        
    

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
    
    nrand=nrng.default_rng(seed)
    
    qList=np.zeros((orderMax,5))
    mList=np.zeros((orderMax*2+2,4))
    mList[1,0]=.00001
    tau,i=FPC.changeTauRe(0,tauMax,mList,pExt,0,mu,m)
    tauList=[tau]
    mList[0:2,0]=[0,tau]
    mList[0,1:4]=[0,0,pExt]
    
    histList=np.zeros((bins,2))
    histList[:,0]=np.linspace(0, tauMax,bins,endpoint=False)
    histList[:,0]+=histList[1,0]*.5
    
    orderList=np.zeros(orderMax+1)
    
    deltaTau=tauMax/bins
    
    total=sum(P)
    pTau=P[0]/total
    pIns=sum(P[:2])/total
    pRem=sum(P[:3])/total
    pSwap=sum(P[:4])/total
    pEx=sum(P[:5])/total
    pFex=sum(P[:6])/total
    #print(pTau,pIns,pRem,pSwap)
    countT=0
    countI=0
    countR=0
    countS=0
    countE=0
    countFE=0
    
    countTD=1
    countID=1
    countRD=1
    countSD=1
    countED=1
    countFED=1
    
    countZero=0
    
    countTherm=1
    countLoopNum=1
    
    
    mcTime=[0,1]
    mcT=1
    n=0
    #startTime=time.time()
    #endTime=runTime*3600+time.time()
    
    
    #change this to a for loop
    for time in range(runTime):
    #while time.time()<endTime:
        #if time.time()-startTime
        x=nrand.uniform()
        #print(pTau<x<pIns,x)
        #print('q',qList)
        #print('m',mList)
        #print(n)
        
        
        #print(x)
        if 0<=x<pTau and n==0:
            #change time zero order
            #print('tau')
            tau,i = FPC.changeTau(tau,tauMax,mList,pExt,n,mu,m)
            
            
            if countLoopNum==step and debug==1:
                
                tauList.append(tau)
                
            countT += i
            countTD+=1
            
            mList[2*n+1,0]=tau
            #print('tau')
            
            #orderList.append(n)  
        elif pTau<x<=pIns and n<orderMax:
            #insert update
            #print('ins')
            
            qList,mList,i=FPC.insertArc(qList,mList,tau,omega,m,n,pIns,pRem,alpha,mu)
            
            countI+=i
            
            countID+=1
            n+=i
            #orderList.append(n)  
            #print(i,n,'in')
            
        elif pIns<x<=pRem and n>=1:
            #remove update
            #print('rem')
            qList,mList,i=FPC.removeArc(qList,mList,omega,m,n,mu,pRem,pIns,alpha)
            
            countR+=i
            
            countRD+=1
            n+=i
            #orderList.append(n)  
            #print(i,n,'rem')
            
        elif pRem<=x<pSwap and n>=2:
            #swap update
            #print('swap')
            qList,mList,i=FPC.swap(qList,mList, n,omega,mu,m)
            countS+=i
            countSD+=1
            #orderList.append(n)  
            
        elif pSwap<=x<pEx and n<=1:
            #extend 

            tau,i = FPC.changeTau(tau,tauMax,mList,pExt,n,mu,m)

            
            if debug==1 and countLoopNum==step:
                tauList.append(tau)
                
                
            countE += i
            countED+=1
            mList[2*n+1,0]=tau
            #print('tau')
            
            #orderList.append(n)  
            
        elif pEx<=x<pFex and n>=1:
            #update is broken and i am too lazy to fix
            #this is a different extend where it rescales the time values relitive to the new tau
            
            qList,mList,tau,i=FPC.fancyExtend(tau,tauMax,mList,qList,pExt,n,mu,m)
            countFE+=i
            countFED+=1
            
            if debug==1 and countLoopNum==step:
                tauList.append(tau)
                
        
        
        
        mcT+=1
        if n==0:
            
            mcTime[0]=(mcTime[0]*mcTime[1]+mcT)/(mcTime[1]+1)
            mcTime[1]+=1
            mcT=0
        
        if debug==1:
            print(n,tau,time)
    
        if thermal<=countTherm and countLoopNum==step:
            #
            #print(n)

            
            if n==0:
                countZero+=1#np.exp(tau*mu)

            histList[int(tau/(deltaTau)),1]+=1#np.exp(tau*mu)

            countLoopNum=0
            
            
            orderList[n]+=1
            
            
        if thermal>countTherm:
            
            countTherm+=1
        else:
            countLoopNum+=1
            

        
        if mcT>mcTMax and mcTMax!=-1:
            #if mcT is set to -1 then this will never happen
            
            qList=np.zeros((orderMax,5))
            mList=np.zeros((orderMax*2+2,4))
            
            tau=FPC.changeTau(0,tauMax,mList,pExt,0,mu,m)
            tauList=[tau[0]]
            mList[0:2,0]=[0,tau[0]]
            mList[0,0:]=pExt
            mcTime.append(mcT)
            countTherm=0
            n=0
            mcT=1
            print('reset')
            
        
    #mcTime.append(mcT)
    count =np.array([mcTime[0],mcTime[1],countT,countI,countID,-countR,countRD,countS,countSD,countE,countED,countFE,countFED])
    if debug==1:
        
        return tauList,countZero,histList,qList,count,mList,orderList
    else:
        
        return histList,countZero,count,orderList

#check out end of second talk
#check how often a update is being rejected 
#check to see if first order is working 

def firstOrderSolution(tau,mu,alpha,omega=1,m=1):
    return -2*alpha*2**.5*np.pi*np.exp(tau*(mu-omega))*m**.5*(2*(omega*tau)**.5+np.exp(omega*tau)*np.pi**.5*(2*omega*tau-1)*erf((omega*tau)**.5))/(32**.5*omega**(1.5)*np.pi**1.5)



    
    
#should I make a plot individual?

def calc(histdata,tauMax,deltaTau,pExt,mu,zeroOrder,m=1,omega=1):
    


    epsK=pExt**2/(2*m)
    
    integral=1/(epsK-mu)*(np.exp(-(epsK-mu)*tauMax)-1)
    

    return -histdata*integral/(deltaTau*(-zeroOrder))
        

        

        


#if i limit the order then i can use numpy array
#this would require the order function to be changed would need to store 
#the order as a seperate variable or can remake the list every time change. 
#can make paralell by using multipul seeds then running the simulation in parallel it looks like


#modulo function
#take floor tau/bin width put in bin 


#data=first_order(20,.001,[100,0,0,0,0,0],0,-6,5,0,1,1,42)

        
        
