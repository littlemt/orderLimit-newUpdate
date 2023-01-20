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
#G_0^tilde=np.exp(-tau*(epsilon-mu))
#epsilon =k**2/(2*m)
#D^tilde =np.exp(-omega*tau)   

def norm(array):
    return np.linalg.norm(array)

def normVec(array):
    dumArray=np.zeros(len(array))
    for i in range(len(array)):
        dumArray[i]=np.linalg.norm(array[i])
    return dumArray

#following are the change in tau Components

def tauReweight(tau,eps_p,mu):
    '''
    

    Parameters
    ----------
    tau : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
        
        
    notes 
    distributions need to test
    power law b*t^c (prob not a polynomial)
    fast growing b*t^(c*t)
    exponential3 b*e^(c*t)
    
    this is probobly good to test first
    these constants should probobly depend on eps_p-mu
    or maybe alpha+mu
    
    do i want it to get just under tau max or do i care?
    
    
    
    e^2*tau does not let it get small enough 
    maybe its rejecting everything that is small?
    scaled(tauNew)/tauOld 
    unless the acceptance ratio is scaled(tauNew)/tauNew
    
    tau^2 returned way to sharp a graph with scaled(tauNew)/tauNew
    
    linear scaling seems to work
    
    need to check acceptance ratio e^(eps-mu)(tauOld-tauNew)
    
    

    '''
    return np.exp(tau*(eps_p-mu))

def tauProb(R,eps_p,mu):
    #need to define the rabdin tau generating function
    return 20*R
    
def changeTau(tau,tauMax,mList,pExt,order,mu,m):
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
    
    eps=pExt**2/(2*m)
    
    
       
        
    t=mList[2*order,0]
    
    
    R=nrand.uniform()
    tauNew=t-np.log(R)/abs(eps-mu)
    
        
        
    
    if tauNew>tauMax or nrand.uniform()>1:
        #np.exp((eps-mu)*(tau-tauNew))*tauScale(tauNew)/tauScale(tau)
        
        return tau,0
    
    else:
        
        return tauNew,1
def changeTauRe(tau,tauMax,mList,pExt,order,mu,m):
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
    
    eps=pExt**2/(2*m)
    R=nrand.uniform()
    
    if order!=0:
        #currently not working need to fix
        
        t=mList[2*order,0]
        
        
        R=nrand.uniform()
        tauNew=t+tauProb(R,eps,mu)
        
        
    else:
        
        tauNew=tauProb(R,eps,mu)
        
    
    if tauNew>tauMax or nrand.uniform()>1:
        #np.exp((eps-mu)*(tau-tauNew))*tauScale(tauNew)/tauScale(tau)
        
        return tau,0
    
    else:
        
        return tauNew,1
    
def fancyExtend(tau,tauMax,mList,qList,pExt,order,mu,m):
    t=mList[2*order,0]
    
    eps=pExt**2/(2*m)
    
    R=nrand.uniform()
    tauNew=t-np.log(R)/abs(eps-mu)
    ratio=tauNew/t
    
    qList[:,0:2]=qList[:,0:2]*ratio
    mList[:2*order+1,0]=mList[:2*order+1,0]*ratio
    
    if tauNew>tauMax:
        return qList,mList,t,0
    else:
        return qList,mList,tauNew,1

#following are functions having to do with insert
    
def spliceInsert(index1,insList,recList,index2):
    length=len(insList)
    
    
    recList[index1+length-2:index2+length-2]=recList[index1:index2]
    
    recList[index1:index1+length]=insList
    return(recList)

def spliceInsertM(index1,insList,recList,index2):
    length=len(insList)
    
    
    recList[index1+length-1:index2+length-1]=recList[index1:index2]
    
    
    recList[index1:index1+length]=insList
    
    return(recList)

    


def insertArc(qList,mList,tMax,orderMax,omega,m,n,pIn,pRem,alpha,mu):
    

    
    index1=nrand.integers(0,2*n+1)
    
    
    #takes the end points of the propogator that was picked at random
    #and gets the endpoints and momentum.
    tauOne=mList[index1,0]
    tauOneP=mList[index1+1,0]
    k=mList[index1,1:]
    
    if tauOneP<tauOne:
        dummy=tauOne
        tauOne=tauOneP
        tauOneP=dummy
        
    #does the calculation of the new arc
    
    #print(mList)
    #print(tauOne,tauOneP)
    tauTwo=nrand.uniform(tauOne,tauOneP)
    
    tauTwoP=tauTwo-np.log(nrand.uniform())/(norm(mList[0,1:])**2/(2*m)-mu)
    
    #checks to see if it is inserting past maxium allowed value
    #print(tauTwoP,tauOneP)
    if tauTwoP>tauOneP:
        return qList,mList,0
    
    
    
    sigma=(m/(tauTwoP-tauTwo))**.5
    #does this have to be abs?
    qTwo=nrand.normal(0,sigma,3)
    
    
    
    
    
    tauListP=np.array([tauOne,tauTwo,tauTwoP,tauOneP])
    momListP=np.ndarray((3,3))
    
    momListP=k,k-qTwo,k
    
    r,this=R_insert(tauListP,momListP,np.array([tauOne,tauOneP]),k,alpha,m,mu,omega,qTwo,pRem,pIn,n)
    
    #print(index1,'i1')
    if r==1:
        #print('i',n)
        qList[n,:2]=[tauTwo,tauTwoP]
        qList[n,2:]=qTwo
        #print(mList)
        
        mList[:,0]=spliceInsert(index1, tauListP, mList[:,0], 2*n+2)
        #print(momListP,index1)
        mList[:,1:]=spliceInsertM(index1, momListP, mList[:,1:], 2*n+1)
        #print(qList,'qL')
        #print(mList,'mL')
        #print('ins')
        
        return qList,mList,1
    else:
        #print('fail')
        return qList,mList,0
        
def R_insert(tauListIn,momentumListIn,tauListRem,momentumListRem,alpha,m,mu,omega,q,pRem,pIn,order):
    #these may be missing |V^2|
    dum1=len(tauListIn)
    #print(tauListIn,'tLI')
    #print('rIn')
    deltaTauListIn=tauListIn[1:dum1]-tauListIn[:dum1-1]
    deltaTauIn=tauListIn[-1]-tauListIn[0]
    #print(deltaTauListIn,deltaTauIn,'dtLI,dtI')
    dum2=len(tauListRem)
    
    
    #print(tauListRem[1:dum2+1],tauListRem[:dum2-1])
    deltaTauListRem=tauListRem[1:dum2+1]-tauListRem[:dum2-1]
    deltaTauRem=tauListIn[-2]-tauListIn[1]
    #print(deltaTauListRem,'dtLR')
    #print(momentumListIn,momentumListRem,'mLI,mLR')
    
    
    alphaTildaSq=2*np.pi*alpha*2**.5
    
    
    wIns=alphaTildaSq*np.exp(-np.sum(deltaTauListIn*(normVec(momentumListIn)**2/2/m-mu)))*np.exp(-omega*(deltaTauIn))*norm(q)**-2*(2*np.pi)**-3
    wRem=np.exp(-np.sum(deltaTauListRem*(normVec(momentumListRem)**2/(2*m)-mu)))

    pXY=pRem*(1/(order+1))
    pYX=pIn/deltaTauIn*omega*np.exp(-omega*(deltaTauRem))*np.exp(-(norm(q)**2/(2*m)*deltaTauIn))\
        /(2*np.pi*m/(deltaTauIn))**(3/2)
    #print(wIns*wRem/pYX*pXY)
    
    if wIns*pYX>wRem*pXY:
        R=1
    elif wIns*pYX>wRem*pXY*1e10:
        R=0
    else:
        R=wIns*pYX/wRem/pXY
    
    
    if nrand.uniform()<R:
        return 1,deltaTauListRem
    else:
        return 0,deltaTauListRem
    
#the following are the functions used in remove
def spliceRemove(index1,remList,recList,index2):
    #print(recList)
    #print(remList,'rL')
    length=len(remList)
    #print(np.shape(remList))
    #print(recList[index1:index1+length],remList)
    
    #this appends the new list to the correct place
    recList[index1:index1+length]=remList
    #print(recList)
    #this takes the end part of the list 
    recList[index1+length:index2-2]=recList[index1+length+2:index2]
    recList[index2-2:]=0
    #print(recList)
    return(recList)
    
def removeArc(qList,mList,orderMax,omega,m,n,mu,pRem,pIn,alpha):
    #print("rem")
    
    #n=order(qList)
    
    #pick random arc
    i=nrand.integers(0,n)
    
    #do i need to check that this does not cross an arc for this update?
    #is this taken into account in the calculation of R
    #the paper says i need to but how does it compare to the video version
    
    
    #the way i am going to do this for now is use the momentum from the bare propogator 
    #that the first vertex is in.

    
    
    
    [tauOne,tauOneP]=qList[i,0:2]
    q=qList[i,2:5]
    #print(qList,n)
    #print(tauOne,tauOneP)
    #print(mList,'m')
    #print(np.where(mList[:,0]==tauOne))
    #print(np.where(mList[:,0]==tauOneP))
    
    #find where the endpoints of the arc are in the mList this gives their time ordered positions
    
    index1=np.where(mList[:,0]==tauOne)[0][0]
    index2=np.where(mList[:,0]==tauOneP)[0][0]
    
    if index2<index1:
        dum=index1
        index1=index2
        index2=dum
    
        
    #if i allow it to pick any arc then I have to make R a large calc
    
    #print(tauOne,tauOneP)    
    r,tauListRem,mListRem =R_remove(qList,mList, index1, index2, m, mu, norm(q), omega, pRem, pIn, n, alpha)
    #print(qList,index1,index2)
    #print(mList,'here')
    if r==1:
        #print('r',n)
        #print(qList)
        #print(mListRem)
        qList[i:n-1]=qList[i+1:n]
        qList[n-1]=np.zeros(5)
        #print(qList)
        #print(mList[:,1])
        mList[:,0]=spliceRemove(index1-1,tauListRem,mList[:,0],2*n+2)
        mList[:,1:]=spliceRemove(index1-1, mListRem, mList[:,1:4], 2*n+1)
        #print(mList)
        #print('rem')
        return qList,mList,-1
    else:
        #print('r')
        return qList,mList,0
    
def R_remove(qList,mList,index1,index2,m,mu,q,omega,pRem,pIn,order,alpha):
    #index1 is the first vertex point
    #index2 is the final xertex 
    #print('rRem')
    momentumList=mList[:,1:4]
    tauList=mList[:,0]
    
    #print(index1,index2)
    
    momentumListIn=momentumList[index1-1:index2+1]
    
    #issue here
    deltaTauListIn=tauList[index1:index2+2]-tauList[index1-1:index2+1]
    deltaTauIn=tauList[index2]-tauList[index1]
    
    deltaTauRem=tauList[index2]-tauList[index1]
    
    #print(tauList)
    dummy=tauList[index1-1:index2+2]
    dummy=np.delete(dummy, (1,-2))
    deltaTauListRem=dummy[1:]-dummy[:-1]
    
    if len(momentumListIn)==3:
        #print('testing')
        momentumListRem=np.array([momentumListIn[0]])
        #print(momentumListRem,1)
        #print(momentumListIn)
        momentumListRemN=norm(momentumListRem)
    else:
        #breaking here when swap is used
        #print(momentumListIn)
        momentumListRem=np.delete(momentumListIn, (1,-2),0)
        momentumListRem[1:-1]+=q
        #print(momentumListRem,2)
        momentumListRemN=normVec(momentumListRem)
    
    alphaTildaSq=2*np.pi*alpha*2**.5
    
    #print(np.shape(deltaTauListIn),np.shape(normVec(momentumListIn)**2/2/m-mu),np.shape(np.exp(-omega*(deltaTauIn))))
    wIns=alphaTildaSq*np.exp(-np.sum(deltaTauListIn*(normVec(momentumListIn)**2/2/m-mu)))*np.exp(-omega*(deltaTauIn))*q**-2*(2*np.pi)**-3
    wRem=np.exp(-np.sum(deltaTauListRem*((momentumListRemN)**2/(2*m)-mu)))

    pXY=pRem*(1/(order+1))
    pYX=pIn/deltaTauIn*omega*np.exp(-omega*(deltaTauRem))*np.exp(-(q**2/(2*m)*deltaTauIn))\
        /(2*np.pi*m/(deltaTauIn))**(3/2)
    
    
    #this should take the tauList and reshape it into an (x,2) array
    tauListRem=tauList[index1-1:index2+2]
    #print(tauListRem)
    tauListRem=np.delete(tauListRem, (1,-2))
    '''tauInd=int(len(tauListRem)/2)
    tauListRem=np.reshape(tauListRem,(tauInd,2))'''
    
    
    
    #print(tauListRem)
    
    #print(momentumListRem,'mR')
    
    #print(wRem/wIns*pYX/pXY)
    
    if wIns*pYX<wRem*pXY:
        R=1
    elif wIns*pXY<wRem*pYX*1e10:
        R=0
    else:
        R=wRem/wIns*pYX/pXY
    
    if nrand.uniform()<R:
        return 1,tauListRem,momentumListRem
    else:
        return 0,tauListRem,momentumListRem
          
        
#the following all has to do with swap   
        
def findEndPoint(qList,tau):
    a=np.where(tau==qList)
    
    if a[1]==0:
        b=a[0],1
    else:
        b=a[0],0
    
    tauP=qList[b]
    q=qList[a[0],2:]
    
    return tauP,q,a,b
    
        
def swap (qList,mList,order,omega,mu,m):
    #updated has bugs
    
    #print('swap')
    
    
    
    #dont remember if integers is inclusive
    a=nrand.integers(2,2*order+1)
    
    tauOne=mList[a,0]
    
    
    #this just picks the closest vertex
    #for this to work it is important that tauOne<tauB
    if a==1:
        tauTwo=mList[2,0]
        b=2
    elif abs(tauOne-mList[a-1,0])>abs(tauOne-mList[a+1,0]):
        
        tauTwo=mList[a+1,0]
        b=a+1
    else:
        
        tauTwo=tauOne
        tauOne=mList[a-1,0]
        b=a
        a=b-1
        
    
    
    tauA,q1,i1,i1p=findEndPoint(qList, tauOne)
    tauB,q2,i2,i2p=findEndPoint(qList, tauTwo)
    
    
    k1=mList[a,1:]
    k1P=swapDecTree(tauOne, tauTwo, tauA, tauB, k1, q1, q2)
    
    x=nrand.uniform()
    
    wX=np.exp(-omega*(abs(tauOne-tauA)+abs(tauTwo-tauB))-(tauTwo-tauOne)*(norm(k1)**2/(2*m)-mu))
    wY=np.exp(-omega*(abs(tauOne-tauB)+abs(tauTwo-tauA))-(tauTwo-tauOne)*(norm(k1P)**2/(2*m)-mu))
    
    if wX>wY:
        r=1
    elif wX>wY*1e10:
        r=0
    else:
        r=wX/wY
    
    if x<r:
        
        qList[i1]=tauTwo
        qList[i2]=tauOne
        
        
        #mList[:,0] remains unchanged
        
        mList[a,1:]=k1P
        
        
        
        return qList,mList,1
    else:
        return qList,mList,0
        
def swapDecTree(t1,t2,ta,tb,k1,q1,q2):
    
    #this function should make the decision as too what type of system it is and outup the correct momentum for the new system
    
    if t2<ta:
        if tb<t1:
            k=k1+q1+q2
        else:
            if ta<tb:
                k=k1+q1
            else:
                k=k1-q1
    else:
        if t1<tb:
            k=k1-q1-q2
        else:
            if ta<tb:
                k=k1+q2
            else:
                k=k1-q2
    return k
            
            

    
    
def changeP(pList):
    #import list of allowed p values
    
    i=nrand.integers(0,len(pList))
    return pList[i]
    
    
    
    

          
    


    
    
    
    
    

    
            
    

#data=first_order(20,1,[100,10,10,.1,1,0],0,-6,5,2,500000)

        