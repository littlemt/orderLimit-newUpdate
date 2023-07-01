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



def normVec(array):
    try:
        return np.linalg.norm(array,axis=1)
    except:
        return np.linalg.norm(array)
   

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
    tauNew=t-np.log(R)*(abs(eps-mu))**-1
    
        

    #the weight may still be 1 need to show
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
    



def insertArc(qList,mList,tMax,omega,m,n,pIn,pRem,alpha,mu):
    '''
    This function is responsible for inserting phonon arcs

    Parameters
    ----------
    qList : numpy array 
        dimension order max by 5.
    mList : numpy array
        dimension 2*order max+2 by 4.
    tMax : TYPE
        DESCRIPTION.
    omega : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.
    pIn : TYPE
        DESCRIPTION.
    pRem : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.
    mu : TYPE
        DESCRIPTION.

    Returns
    -------
    qList : TYPE
        DESCRIPTION.
    mList : TYPE
        DESCRIPTION.
    int
        DESCRIPTION.

    '''

    
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
    
    
    
    #should be log/omega not this
    tauTwoP=tauTwo-np.log(nrand.uniform())/omega#(norm(mList[0,1:])**2/(2*m)-mu)
    
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
    
    
    # wIns=alphaTildaSq*np.exp(-np.sum(deltaTauListIn*(normVec(momentumListIn)**2/2/m-mu)))*np.exp(-omega*(deltaTauRem))*normVec(q)**-2*(2*np.pi)**-3
    # wRem=np.exp(-np.sum(deltaTauListRem*(normVec(momentumListRem)**2/(2*m)-mu)))
    
    

    # pYX=pRem*(1/(order+1))
    # pXY=pIn/deltaTauIn*omega*np.exp(-omega*(deltaTauRem))*np.exp(-(normVec(q)**2/(2*m)*deltaTauRem))\
    #     /(2*np.pi*m/(deltaTauRem))**(3/2)
        
    wExp=-np.sum(deltaTauListIn*(normVec(momentumListIn)**2/2/m-mu))+np.sum(deltaTauListRem*(normVec(momentumListRem)**2/(2*m)-mu))\
        +(normVec(q)**2/(2*m)*deltaTauRem)
            
    wRatio=pRem*(1/(order+1))*alphaTildaSq*normVec(q)**-2*(2*np.pi)**-3/(pIn/deltaTauIn*omega/(2*np.pi*m/(deltaTauRem))**(3/2))
    
    
    if (wExp)>np.log(1/wRatio):
        R=1
        
    elif wExp<np.log(1E-16/wRatio):
        R=0
        
    else:
        R=(wRatio*np.exp(wExp))
    
    #print(R,wIns*pYX/(wRem*pXY),'i')
    
    if nrand.uniform()<R:
        
        return 1,deltaTauListRem
    else:
        return 0,deltaTauListRem
    
def spliceInsert(index1,insList,recList,index2):
    '''
    function inserts the new verticies and all changed propogators in the appropriate
    time ordered position in the list

    Parameters
    ----------
    index1 : TYPE
        DESCRIPTION.
    insList : TYPE
        DESCRIPTION.
    recList : TYPE
        DESCRIPTION.
    index2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    length=len(insList)
    
    
    recList[index1+length-2:index2+length-2]=recList[index1:index2]
    
    recList[index1:index1+length]=insList
    return(recList)

def spliceInsertM(index1,insList,recList,index2):
    '''
    

    Parameters
    ----------
    index1 : TYPE
        DESCRIPTION.
    insList : TYPE
        DESCRIPTION.
    recList : TYPE
        DESCRIPTION.
    index2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    length=len(insList)
    
    
    recList[index1+length-1:index2+length-1]=recList[index1:index2]
    
    
    recList[index1:index1+length]=insList
    
    return(recList)

    
    
#the following are the functions used in remove

def removeArc(qList,mList,omega,m,n,mu,pRem,pIn,alpha):
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
    r,tauListRem,mListRem =R_remove(qList,mList, index1, index2, m, mu, normVec(q), omega, pRem, pIn, n, alpha)
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
    
    #seperates the electron array into just verticies and just momentums
    momentumList=mList[:,1:4]
    tauList=mList[:,0]
    
   
    
    #takes the two indicies given to it by the remove arc function and creates the array of 
    #the momentums of all of the electron propogators covered by the sugested arc to be removed
    #as well as the electron propogators to the immidiate either side
    momentumListIn=momentumList[index1-1:index2+1]
    
    #this calculates the change in time for each electron propagator
    deltaTauListIn=tauList[index1:index2+2]-tauList[index1-1:index2+1]
    
    
    
    #calculates delta time for the 
    deltaTauIn=tauList[index2+1]-tauList[index1-1]
    
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
        momentumListRemN=normVec(momentumListRem)
    else:
        #breaking here when swap is used
        #print(momentumListIn)
        momentumListRem=np.delete(momentumListIn, (1,-2),0)
        momentumListRem[1:-1]+=q
        #print(momentumListRem,2)
        momentumListRemN=normVec(momentumListRem)
    
    alphaTildaSq=2*np.pi*alpha*2**.5
    
    #print(np.shape(deltaTauListIn),np.shape(normVec(momentumListIn)**2/2/m-mu),np.shape(np.exp(-omega*(deltaTauIn))))
    wIns=alphaTildaSq*np.exp(-np.sum(deltaTauListIn*(normVec(momentumListIn)**2/2/m-mu)))*np.exp(-omega*(deltaTauRem))*q**-2*(2*np.pi)**-3
    wRem=np.exp(-np.sum(deltaTauListRem*((momentumListRemN)**2/(2*m)-mu)))
    
    pYX=pRem*(1/(order))
    pXY=pIn/deltaTauIn*omega*np.exp(-omega*(deltaTauRem))*np.exp(-(q**2/(2*m)*deltaTauRem))/(2*np.pi*m/(deltaTauRem))**(3/2)
        
   
    wExp=-np.sum(deltaTauListRem*((momentumListRemN)**2/(2*m)-mu))\
            +(np.sum(deltaTauListIn*(normVec(momentumListIn)**2/2/m-mu))-(q**2/(2*m)*deltaTauRem))
            
    wRatio=pIn/deltaTauIn*omega*(2*np.pi*m/(deltaTauRem))**(-3/2)/(pRem*(1/(order)))/alphaTildaSq*q**2*(2*np.pi)**3
    
    # if wRatio<0:
    #     print(deltaTauIn,deltaTauRem)
    #     print(tauList)
    #     print(index1,index2)
    #     print(tauList[index2+1],tauList[index1-1])
    #     exit()
   
    #print(R,(wRatio*np.exp(wExp)),order)


    
    
    
    #this should take the tauList and reshape it into an (x,2) array
    tauListRem=tauList[index1-1:index2+2]
    #print(tauListRem)
    tauListRem=np.delete(tauListRem, (1,-2))
    '''tauInd=int(len(tauListRem)/2)
    tauListRem=np.reshape(tauListRem,(tauInd,2))'''
    
    
    
    #print(tauListRem)
    
    #print(momentumListRem,'mR')
    
    
    
    if (wExp)>np.log(1/wRatio):
        R=1
        
    elif wExp<np.log(1E-16/wRatio):
        # print(wExp,np.log(1E-16/wRatio),wRatio)
        R=0
        
    else:
        R=wRatio*np.exp(wExp)
        
        
    # print(R,(wIns*pYX/(wRem*pXY))**-1,'r',wRatio)
   
    if nrand.uniform()<R:
    
       
        return 1,tauListRem,momentumListRem
    else:
        return 0,tauListRem,momentumListRem
          
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
    
        
#the following all has to do with swap   
        
def findEndPoint(qList,tau):
    '''
    this functuion is used to find the appropriate phonon arc for the chosen end vetex
    

    Parameters
    ----------
    qList : numpy array
        array contaning all the phonon arcs.
    tau : float
        chosen vertex.

    Returns
    -------
    tauP : float
        the other vertex in the phonon arc.
    q : numpy array
        the momentum of the phonon arc.
    a : float
        the index of the chosen vertex.
    b : float
        the index of the other vertex.

    '''
    #print(tau,qList,'start')
    a=np.where(tau==qList)
    
    
    #this finds the other end point since it can only be in positions 0 or 1 it checks and assigns the other

    if a[1]==0:
        #print('a1')
        b=a[0],1
    else:
        #print('b2')
        b=a[0],0
    
    
    tauP=qList[b]
    q=qList[a[0],2:]
    
    return tauP,q,a,b
    
        
def swap (qList,mList,order,omega,mu,m):

    
    #print('swap')
    
    
    
    a=nrand.integers(1,2*order+1)
    
    tauOne=mList[a,0]
    
    
    #this just picks the closest vertex 
    #for this to work it is important that tauOne<tauTwo
    if a==1:
        tauTwo=mList[2,0]
        b=2
    elif a==2*order:
        tauTwo=mList[2*order,0]
        tauOne=mList[2*order-1,0]
        a=2*order-1
        b=2*order
    elif abs(tauOne-mList[a-1,0])>abs(tauOne-mList[a+1,0]):
        
        tauTwo=mList[a+1,0]
        b=a+1

    else:
        
        tauTwo=tauOne
        tauOne=mList[a-1,0]
        b=a
        a=b-1
        
    
    #this section finds the other vetex for both of the chosen vertex
    tauA,q1,i1,i1p=findEndPoint(qList, tauOne)
    tauB,q2,i2,i2p=findEndPoint(qList, tauTwo)
    
    
    if tauTwo==tauA or tauTwo==tauOne:
        return qList,mList,0
    #print(findEndPoint(qList, tauTwo))
    k1=mList[a,1:]
    k1P=swapDecTree(tauOne, tauTwo, tauA, tauB, k1, q1, q2)
    
    x=nrand.uniform()
    
    
    
    
    xExp=-omega*(abs(tauOne-tauA)+abs(tauTwo-tauB))-(tauTwo-tauOne)*(normVec(k1)**2/(2*m)-mu)
    yExp=-omega*(abs(tauOne-tauB)+abs(tauTwo-tauA))-(tauTwo-tauOne)*(normVec(k1P)**2/(2*m)-mu)
    # wX=np.exp(-omega*abs(tauOne-tauA))*np.exp(-omega*abs(tauTwo-tauB))*np.exp(-(tauTwo-tauOne)*(normVec(k1)**2/(2*m)-mu))
    # wY=np.exp(-omega*abs(tauOne-tauB))*np.exp(-omega*abs(tauTwo-tauA))*np.exp(-(tauTwo-tauOne)*(normVec(k1P)**2/(2*m)-mu))
    
    
    
    
    # print(tauOne,tauTwo,tauA,tauB)
    # print(qList[i1[0]],'arc1')
    # print(qList[i2[0]],'arc2')
    # print(k1,k1P,'m')
    # print(np.exp(xExp-yExp),wX/wY,'exponents')
    # if wX==0 or wY==0:
    #     print(wX,wY,'wx/y')
    #     print(k1P)
    #     print(qList,mList,'arrays')
    #     print(tauOne,tauA,tauTwo,tauB,q1,q2,'index')
    #     print(np.exp(xExp-yExp))
    #     exit()
        
    
    if xExp-yExp>0:
        # print(np.exp(xExp-yExp),1,'r')
        r=1
    elif xExp-yExp<-36:
        r=0
        # print(np.exp(xExp-yExp),r,'r')
    else:
        r=np.exp(xExp-yExp)
        # print(np.exp(xExp-yExp),r,'r')
    
    
    if x<r:
        
        #takes the index given from find and swaps the positions of the two 
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
            #print(1)
            #case 6 
            k=k1+q1+q2
        else:
            if ta<tb:
                #case 1
                #print(2)
                k=k1+q1-q2
            else:
                #case 2
                #print(3)
                k=k1+q1-q2
    else:
        if t1<tb:
            #print(4)
            #case 5
            k=k1-q1-q2
        else:
            if ta<tb:
                #case 4
                #print(5)
                k=k1+q2-q1
            else:
                #case 3
                #print(6)
                k=k1-q1+q2
    return k
            
            

    
    
def changeP(pList):
    #import list of allowed p values
    
    i=nrand.integers(0,len(pList))
    return pList[i]
    
    
    
    

          
    


    
    
    
    
    

    
            
    

#data=first_order(20,1,[100,10,10,.1,1,0],0,-6,5,2,500000)

        