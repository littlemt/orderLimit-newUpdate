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
    
    

def insertArc(qList,tMax,orderMax,omega,m,n,pIn,pRem,alpha,mu):
    
    print("ins")
    
    
    index1=nrand.integers(0,n+1)
    if index1==n:
        index2=0
    else:
        index2=nrand.integers(0,2)
    
    
    dummyList=np.zeros((n+1-index1,4))
    tauList=qList[index1:n+1,3:5].flatten()
    momList=qList[index1:n+1,5:7].flatten()
    
    
    if index2==0:
        print("i0")
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
        #fix this so -ln(x)/scaling
        tauTwoP=tauTwo+nrand.exponential(1/omega,None)
        
        #checks to see if it is inserting past maxium allowed value
        
        if tauTwoP>tMax:
            return qList,0
        
        dTa=tauTwoP-tauTwo
        
        sigma=(m/(tauTwoP-tauTwo))**.5
        qTwo=nrand.normal(0,sigma)
        
        #this makes a dummy list that contains the bare propogators in time order
        #inserts the new values then moves the other values down the list 
        
        
        
        for i in range(index1*2+index2,len(tauList)):
            if tauList[i]<tauTwoP<tauList[i+1]:
                dumInd=i
                break
            
        
        index1P=[np.floor(dumInd/2),dumInd%2]
        
        dummyList[0]=[tauOne,tauTwo,k,k-qTwo]
        
        momListP=np.array([k,momList[index1*2+index2:dumInd+2]-qTwo,momList[dumInd+1],0]).flatten()
        tauListP=np.array([tauOne,tauTwo,tauList[index1*2+index2+1:dumInd+1],tauTwoP,tauOneP]).flatten()
        
        
        dummyList[1:index1P[0]+1-index1[0],2:4]=np.reshape(momListP,(dumInd+1-index1,2))
        dummyList[1:index1P[0]+1-index1[0],:2]=np.reshape(tauListP,(dumInd+1-index1,2))
        
        #lets hope everything is sliced right
        
        
        
        #all comments apply to each associated section in the else
        
    else:
        print(1)
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
        
        
         
        dumInd=index1*2+index2
        while tauTwoP>tauList[dumInd]:
            dumInd+=1
        
        index1P=[np.floor(dumInd/2),dumInd%2]
        
        dummyList[0]=[tauOne,tauTwo,k,k-qTwo]
        
        momListP=np.array([k,momList[index1*2+index2:dumInd+2]-qTwo,momList[dumInd+1]]).flatten()
        tauListP=np.array([tauOne,tauTwo,tauList[index1*2+index2+1:dumInd+1],tauTwoP,tauOneP]).flatten()
        
        
        dummyList[1:index1P[0]+1-index1[0],2:4]=np.reshape(momListP,(dumInd-index1,2))
        dummyList[1:index1P[0]+1-index1[0],:2]=np.reshape(tauListP,(dumInd-index1,2))
        
    
    
    
    
    
    r=R_insert(tauListP,momListP,tauList[index1*2+index2:dumInd+1],momList[index1*2+index2:dumInd+1],alpha,m,mu,omega,qTwo,pRem,pIn,n)
    
    x=nrand.uniform()
    
    if x<r:
        qList[n,0:3]=np.array([tauTwo,tauTwoP,qTwo])
        qList[index1:n+1,3:7]=dummyList
       
         
        #need to check if i am splicing right
        return qList,1
    else:
        return qList,0
        
    
    
    
def removeArc(qList,omega,tMax,orderMax,m,n,mu,pRem,pIn,alpha):
    print("rem")
    
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
    
    '''if index1[1]==0:
        tauTwo=qList[index1[0]-1,1]
        
        
    else:
        tauTwo=qList[index1[0],0]
        
    if index2[1]==0:
        tauTwoP=qList[index2[0]-1,1]
    else:
        tauTwoP=qList[index2[0],0]'''
        
    #if i allow it to pick any arc then I have to make R a large calc
    
        
    r,propList=R_remove(qList, index1, index2, m, mu, q, omega, pRem, pIn, n, alpha)
    
    if r==1:
        #print('a',i)
        if i==n-1:
            #this part says if the arc to remove is the same as the max order 
            #then just remove the last arc from the list 
            qList[n-1,:3]=np.zeros(3)
            qList[index1[0]-1:index2[0]+1,3:7]=propList
        else:
            #takes the arcs on the list from after the arc picked and places them back one spot
            
            dummy=qList[i+1:orderMax]
            qList[i:orderMax-1]=dummy
            qList[orderMax-1]=np.zeros(3)
            
            qList[index1[0]-1:index2[0]+1,3:7]=propList
        
        return qList,-1
    else:
        #print('r')
        return qList,0
        
        
        
def swap (qList,order):
    #not updated yet
    
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
    
    
    
    
def R_insert(tauListIn,momentumListIn,tauListRem,momentumListRem,alpha,m,mu,omega,q,pRem,pIn,order):
    
    dum1=len(tauListIn)
    deltaTauListIn=tauListIn[1:dum1+1]-tauListIn[:dum1]
    deltaTauIn=tauListIn[-1]-tauListIn[0]
    
    dum2=len(tauListRem)
    deltaTauListRem=tauListRem[1:dum2+1]-tauListRem[:dum1]
    deltaTauRem=tauListIn[-2]-tauListIn[1]
    
    
    alphaTildaSq=2*np.pi*alpha*2**.5
    
    ratio1=alphaTildaSq*np.exp(-np.sum(deltaTauListIn*(momentumListIn/2/m-mu)))*np.exp(omega*(deltaTauIn))\
        / np.exp(-np.sum(deltaTauListRem*((momentumListRem)/(2*m)-mu)))*q**2
        
    ratio2=pRem/pIn*(1/(order))*deltaTauIn*(2*m*np.pi/deltaTauIn)**(3/2)\
        /(omega*np.exp(-omega*(deltaTauRem))*np.exp(-q**2/(2*m)*deltaTauRem))
          
    if nrand.uniform()<1/ratio1/ratio2:
        return 1,deltaTauListRem
    else:
        return 0,deltaTauListRem
          
    

def R_remove(qList,index1,index2,m,mu,q,omega,pRem,pIn,order,alpha):
    
    
    momentumList=qList[:,-2:]
    tauList=qList[:,3:5]
    
    momentumList=momentumList.flatten()
    tauList=tauList.flatten()
    
    
    
    i1=index1[0]*2+index1[1]
    i2=index2[0]*2+index2[1]
    
    momentumListIn=momentumList[i1-1:i2+2]
    deltaTauListIn=tauList[i1:i2+3]-tauList[i1-1:i2+2]
    deltaTauIn=tauList[i2+1]-tauList[i1-1]
    deltaTauRem=tauList[i2]-tauList[i1]
    
    
    dummy=deltaTauListIn[1],deltaTauListIn[-2]
    deltaTauListRem=np.delete(deltaTauListIn, (1,-2))
    deltaTauListRem[0],deltaTauListRem[-1]=deltaTauListRem[0]+dummy[0],deltaTauListRem[-1]+dummy[1]
    
    momentumListRem=np.delete(momentumListIn, (1,-2))
    momentumListRem[1:-1]+=q
    
    alphaTildaSq=2*np.pi*alpha*2**.5
    
    ratio1=alphaTildaSq*np.exp(-np.sum(deltaTauListIn*(momentumListIn/2/m-mu)))*np.exp(omega*(deltaTauIn))\
        / np.exp(-np.sum(deltaTauListRem*((momentumListRem)/(2*m)-mu)))*q**2
        
    ratio2=pRem/pIn*(1/(order))*deltaTauIn*(2*m*np.pi/deltaTauIn)**(3/2)\
        /(omega*np.exp(-omega*(deltaTauRem))*np.exp(-q**2/(2*m)*deltaTauRem))
          
    if nrand.uniform()<1/ratio1/ratio2:
        return 1,deltaTauListRem
    else:
        return 0,deltaTauListRem
          
    
    
    
    
    
    
'''    
def order(qList):
    x=np.count_nonzero(qList[0])
    return x
    '''
    
            
    
#might need to fix update insert remove

        