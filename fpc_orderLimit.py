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


def first_order(tauMax,runTime,P,pExt,mu,alpha,orderMax,thermal,step,mcTMax=-1,bins=100,omega=1,m=1,debug=0):
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
    
    
    qList=np.zeros((orderMax,5))
    mList=np.zeros((orderMax*2+2,4))
    
    tau,i=FPC.changeTau(0,tauMax,mList,pExt,0,mu,m)
    tauList=[tau]
    mList[0:2,0]=[0,tau]
    mList[0,1:4]=[0,0,pExt]
    
    histList=np.zeros((bins,2))
    histList[:,0]=np.linspace(0, tauMax,bins)
    
    deltaTau=tauMax/bins
    
    total=sum(P)
    pTau=P[0]/total
    pIns=sum(P[:2])/total
    pRem=sum(P[:3])/total
    pSwap=sum(P[:4])/total
    pExt=sum(P[:5])/total
    pFext=sum(P[:6])/total
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
    
    #orderList=[]
    mcTime=[0,1]
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
            
            qList,mList,i=FPC.insertArc(qList,mList,tau,orderMax,omega,m,n,pIns,pRem,alpha,mu)
            
            countI+=i
            countID+=1
            n+=i
            #orderList.append(n)  
            #print(i,n,'in')
            
        elif pIns<x<=pRem and n>=1:
            #remove update
            #print('rem')
            qList,mList,i=FPC.removeArc(qList,mList,orderMax,omega,m,n,mu,pRem,pIns,alpha)
            
            countR+=i
            countRD+=1
            n+=i
            #orderList.append(n)  
            #print(i,n,'rem')
            
        elif pRem<=x<pSwap and n>=2:
            #swap update
            qList,mList,i=FPC.swap(qList,mList, n,omega,mu,m)
            countS+=i
            countSD+=1
            #orderList.append(n)  
            
        elif pSwap<=x<pExt and n<=1:
            #extend 
            tau,i = FPC.changeTau(tau,tauMax,mList,pExt,n,mu,m)
            
            if debug==1 and countLoopNum==step:
                tauList.append(tau)
                
                
            countE += i
            countED+=1
            mList[2*n+1,0]=tau
            #print('tau')
            
            #orderList.append(n)  
            
        elif pExt<=x<pFext and n<=1:
            #this is a different extend where it rescales the time values relitive to the new tau
            
            qList,mList,tau,i=FPC.fancyExtend(tau,tauMax,mList,qList,pExt,n,mu,m)
            countFE+=i
            countFE+=1
            
            if debug==1 and countLoopNum==step:
                tauList.append(tau)
                
                
            
            #orderList.append(n)  
            
            
            
        
        
        
        #orderList.append(n)  
        
        if debug==1:
            mcT+=1
            if n==0:
                
                mcTime[0]=(mcTime*mcTime[1]+mcT)/(mcTime+1)
                mcTime[1]+=1
        
        
        if thermal<=countTherm and countLoopNum==step:
            #
            
            if n==0:
                countZero+=1
            
            histList[int(tau/(deltaTau)),1]+=1
            
            countLoopNum=0
            
        if thermal>countTherm:
            
            countTherm+=1
        else:
            countLoopNum+=1
            
            
        
        if mcT>mcTMax and mcTMax!=-1:
            #if mcT is set to -1 then this will never happen
            print('reset')
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
            
        
    #mcTime.append(mcT)
    count =np.array([countT/countTD,countI/countID,-countR/countRD,countS/countSD,countE/countED,countFE/countFED])
    if debug==1:
        return tauList,mcTime,countZero,histList,qList,count,mList
    else:
        return histList,countZero,count

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

def calc(data,pExt,mu,zeroOrder,m=1,omega=1):
    
    histdata,zeroOrder=data
    tauMax=histdata[-1,0]
    
    deltaTau=data[0,1]

    epsK=pExt**2/(2*m)
    
    integral=1/(epsK-mu)*(np.exp(-(epsK-mu)*tauMax)-1)
    
    
        #just realized this is implimented wrong i think
    histdata[:,2]=-histdata[:,1]*integral/deltaTau/-zeroOrder

    return histdata     
        
def saveData(data,path,tauMax,runTime,P,pExt,mu,alpha,orderMax,therm,step,bins,seed,mctMax=-1):
    dumString='tM'+str(tauMax)+'rT'+str(runTime)+'hr'+'prob'+str(P)+'mom'+str(pExt)\
        +'mu'+str(mu)+'a'+str(alpha)+'oM'+str(orderMax)+'therm'+str(therm)+'step'+str(step)+'bins'+str(bins)+'seed'+str(seed)
    if mctMax!=-1:
        dumString+='lim'+mctMax
    np.savetxt(path+'histList'+dumString,data[0])
    
    np.savetxt(path+'countPercent'+dumString,data[2])
    #np.savetxt(path+'orderList'+dumString,data[3])
    return
        
def histogram(data,title,xAxis,yAxis,scale,binNo,sample):
    mpl.title(title)
    mpl.yscale(scale)
    mpl.xlabel(xAxis)
    mpl.ylabel(yAxis)
    mpl.hist(data,bins=binNo)
    mpl.plot()
    return
        

        
#add extend 



#if i limit the order then i can use numpy array
#this would require the order function to be changed would need to store 
#the order as a seperate variable or can remake the list every time change. 
#can make paralell by using multipul seeds then running the simulation in parallel it looks like


#modulo function
#take floor tau/bin width put in bin 


#note average mcTime between order 0

        
        