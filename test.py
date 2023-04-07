#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:19:13 2022

@author: littlemt
"""

import configparser
import numpy as np
import matplotlib.pyplot as mpl
from  math import floor 
from scipy.integrate import nquad 
import dep.FPDMCUpdates_OrderLimit as FPC
from numpy.linalg import norm

#fix extend 
def run(seed):

    hist,zero,count,order=FPC.first_order(15,100000000,[1000,10,10,10,10,0],0,-6,5,500,100000,100,seed)

    hist2=np.zeros((100,3))
    hist2[:,0]=hist[:,0]
    hist2[:,1]=FPC.calc(hist[:,1],hist[-1,0],hist[1,0]*2,0,-6,zero)
    mpl.scatter(hist2[:,0],np.log(-hist2[:,1]))
    
    mpl.show()
    mpl.bar(np.arange(len(order)),order/sum(order))
    mpl.show()
    print(count)
    print(hist)


def G(tau,p,mu):
    return -np.exp(-tau*(p**2/2-mu))

def D(q,tau):
    return np.exp(tau)/q**2

def alphaTilde(alpha):
    return 2*np.pi*alpha*2**.5
    

def rIns(aList,eList,alpha,propArc,mu):
    wX=-G(eList[1,0],norm(eList[0,1:]),mu)
    wY=alphaTilde(alpha)*G(propArc[0],norm(eList[0,1:]),mu)*G(propArc[1]-propArc[0],norm(eList[0,1:]-propArc[2:]),mu)*\
        G(eList[1,0]-propArc[1],norm(eList[0,1:]),mu)*D(norm(propArc[2:]),propArc[1]-propArc[0])*(2*np.pi)**-3
    pXY=1/(eList[1,0])*np.exp(propArc[1]-propArc[0])*G(propArc[1]-propArc[0],norm(propArc[2:]),0)/(2*np.pi/(propArc[1]-propArc[0]))**(3/2)
    pYX=1
    
    print(wX,wY,pXY,pYX)

    
    return wY/wX*pYX/pXY


    
def rRem(aList,eList,alpha,propArc,mu):
    wX=-G(eList[1,0],norm(eList[0,1:]),mu)
    wY=-alphaTilde(alpha)*G(propArc[0],norm(eList[0,1:]),mu)*G(propArc[1]-propArc[0],norm(eList[0,1:]-propArc[2:]),mu)*\
        G(eList[1,0]-propArc[1],norm(eList[0,1:],mu))*D(propArc[2:],propArc[1]-propArc[0])*(2*np.pi)**-2
    pXY=1/(propArc[1]-propArc[0])*np.exp(propArc[1]-propArc[0])*G(propArc[1]-propArc[0],norm(propArc[2:]),0)/(2*np.pi/(propArc[1]-propArc[0]))**(3/2)
    pYX=1
    
    return wX/wY*pXY/pYX

aList=np.zeros((2,5))
eList=np.zeros((6,4),dtype=float)
eList[1,0]=1
n=0
while n<2:
    aList,eList,i=FPC.insertArc(aList, eList, 5, 1, 1, n, 1, 1, 5, -6)
    n+=i
  

for i in range(100000):    
    FPC.swap(aList,eList,n,1,-6,1)












#print(rIns(aList,eList,5,prop[0],-6),'test')


#rIns(aList,eList,5,)
    