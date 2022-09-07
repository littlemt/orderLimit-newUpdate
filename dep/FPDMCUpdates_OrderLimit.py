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
    

def changeTau(tau,tauMax,qList,pExt,order,mu,m):
    '''
    this is a combination of the change tau and extend built into one

    Parameters
    ----------
    tau : TYPE
        DESCRIPTION.
    tauMax : TYPE
        DESCRIPTION.
    qList : TYPE
        DESCRIPTION.
    pExt : TYPE
        DESCRIPTION.
    order : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    int
        DESCRIPTION.

    '''
    
    t=tau
    if order!=0:
        tList=qList[:order,0:2]
        tList=tList.flatten()
        
        t=np.max(tList)
        
        eps=pExt**2/(2*m)
        
        R=nrand.uniform()
        tauNew=t-np.log(R)/abs(eps-mu)
        
        
    else:
        eps=pExt**2/(2*m)
        
        R=nrand.uniform()
        tauNew=-np.log(R)/abs(eps-mu)
    if tauNew>tauMax:
        return tau,0
    else:
        return tauNew,1
    
    

def insertArc(qList,p,tMax,orderMax,omega,m,n):
    
    
    
    index1=nrand.integers(0,n+1)
    if index1==n:
        index2=0
    else:
        index2=nrand.integers(0,2)
    
    
    dummyList=np.zeros((n+1-index1,4))
    
    if index2==0:
        
        #either part of the if statment does the same thing but it depends which bare 
        #propogator is picked based on the data structure. 
        #qList is a terrible name dont want to change it now though.
        
        
        #takes the end points of the propogator that was picked at random
        #and gets the endpoints and momentum.
        tauOne=qList[index1,index2+3]
        tauOneP=qList[index1,4]
        k=qList[index1,5]
        
        
        #does the calculation of the new arc
        tauTwo=nrand.uniform(tauOne,tauOneP)
        tauTwoP=tauTwo+nrand.exponential(1/omega,None)
        
        if tauTwoP>tMax:
            return qList,0
        
        dTa=tauTwoP-tauTwo
        
        sigma=(m/(tauTwoP-tauTwo))**.5
        qTwo=nrand.normal(0,sigma)
        
        #this makes a dummy list that contains the bare propogators in time order
        #inserts the new values then moves the other values down the list 
        dummyList[0:2]=np.array([[tauOne,tauTwo,k,k-qTwo],[tauTwoP,tauOneP,k,qList[index1,6]]])
        dummyList[2:n+2-index1]=qList[index1+1:n+2,3:7]#unsure about n+1 or n+2
        
        #all comments apply to each associated section in the else
        
    else:
        tauOne=qList[index1,index2+3]
        tauOneP=qList[index1+1,3]
        k=qList[index1,6]
        
        tauTwo=nrand.uniform(tauOne,tauOneP)
        tauTwoP=tauTwo+nrand.exponential(1/omega,None)
        
        if tauTwoP>tMax:
            return qList,0
        
        dTa=tauTwoP-tauTwo
        
        sigma=(m/(tauTwoP-tauTwo))**.5
        qTwo=nrand.normal(0,sigma)
        
        
        
        dummyList[0:2]=[qList[n],[tauTwo,tauTwoP,k-qTwo,k]]
        dummyList[2:n+2-index1]=qList[index1+1:n+2,3:7]
        
    
    
    
    
    
    r=R_insert(omega, m, k, tauTwo, tauTwoP, qTwo, n, dTa)
    
    x=nrand.uniform()
    
    if x<r:
        qList[n,0:3]=np.array([tauTwo,tauTwoP,qTwo])
        qList[index1:n+1,3:7]=dummyList
       
         
        #need to check if i am splicing right
        return qList,1
    else:
        return qList,0
        
    
    
    
def removeArc(qList,omega,tMax,orderMax,m,p,n):
    #n=order(qList)
    
    #pick random arc
    i=nrand.integers(0,n)
    
    #do i need to check that this does not cross an arc for this update?
    #is this taken into account in the calculation of R
    #the paper says i need to but how does it compare to the video version
    
    
    #the way i am going to do this for now is use the momentum from the bare propogator 
    #that the first vertex is in.

    
    
    
    [tauOne,tauOneP,q]=qList[i,0:3]
    
    index1=np.where(qList[:,3:5]==tauOne)
    index2=np.where(qList[:,3:5]==tauOneP)
    
    if index1[1]==0:
        tauTwo=qList[index1[0]-1,1]
        
        
    else:
        tauTwo=qList[index1[0],0]
        
    if index2[1]==0:
        tauTwoP=qList[index2[0]-1,1]
    else:
        tauTwoP=qList[index2[0],0]
        
    #if i allow it to pick any arc then I have to make R a large calc
    #R should be
        
    
    
    
    x=nrand.uniform()
    #print(i,n)
    dum1=qList[i+1:n+1]
    dum2=qList[i:n]
    #print(np.shape(dum1),np.shape(dum2))
    if x<r:
        #print('a',i)
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
        #print('r')
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
    
    
def changeP(pList):
    #import list of allowed p values
    
    i=nrand.integers(0,len(pList))
    return pList[i]
    
    
    
    
def R_insert(omega,m,p,tauTwo,tauTwoP,qTwo,n,dTa):
    #D_n+1 /D_n /Omega *(2n+1)deltaTau/
    #Omega=
    #D_n+1/D_n=
    
    '''
    q=np.linalg.norm(q)
    rad=q**2/2/m*(tauTwo-tauOne)
    return 2*(2*np.pi)**.5*alpha**2*np.exp(rad)*prem*(tauL-tauR)\
            /((1+order)*pins*q**2*omegaPH)*(m/(tauTwo-tauOne))**(3/2)
    '''
    
    D=np.exp(-omega*(tauTwoP-tauTwo))*np.exp(-(qTwo**2-2*p*qTwo)/(2*m)*(tauTwoP-tauTwo))
    qNot=np.sqrt(2*m*omega)
    Omega=1/(4*np.pi)/qNot*omega*np.exp(-omega*(1+qTwo/qNot)**2*(tauTwoP-tauTwo))
    return D/Omega*(2*n+1)*(2*n+1)*dTa/(n+1)

def R_remove(qList,tauOne,tauOneP,m,mu):
    
    '''
    (pow(sqrt(mass / (2*M_PI*dt)), dim) * vertex_amplitude_sq
                    * nr_int * tau_intv * update_prob[remove]/update_prob[insert]
                    * exp(dt* (p_ext* q)/mass) / (omega_p * order * q.r_sq()));
    
    
    nr_int=2*n-1
    tau_intv=
    dt = it2->time() - it1->time()
    vertex_amplitude_aq=alpha^2????
    
    ((m/(2*M_PI*dt))^.5)^dim
    what is M_Pi? m*pi?
    
    
    '''
    eps_mu=qList[tauOne:tauOneP+1,-2:]/2/m-mu
    
    
    
    
'''    
def order(qList):
    x=np.count_nonzero(qList[0])
    return x
    '''
    
            
    
#might need to fix update insert remove

        