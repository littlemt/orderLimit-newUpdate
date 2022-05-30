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
    

def changeTau(tMax,p,m):
    eps=p**2/(2*m)
    
    R=nrand.uniform()
    tauNew=tMax-np.log(R)/eps
    
    return tauNew
    
    

def insertArc(qList,tMax,orderMax):
    
    tauList=qList[0:1]
    tauList=np.trim_zeros(tauList)
    tauList=np.sort(tauList,axis=None)
    
    index=nrand.integer(0,tauList.size)
    
    tauOne=tauList[index]
    tauOneP=tauList[index+1]
    tauTwo=nrand.uniform(tauOne,tauOneP)
    
    #What is omega?
    #distrobution is expotneital or power law dont know how to scale it.  maybe omega?
    
    #p=1/(2*order+1)/(deltaTau)
    
def removeArc(qList,order):
    
    nrand.integers(0,order)
    
    
    
def R_insert():
    #D_n+1 /D_n /Omega *(2n+1)deltaTau/
    #Omega=
    #D_n+1/D_n=
    
    
            
    
    #need to create list of nodes
    
    #new change tau same as old but can do at not 0 order 1:03:40
    #new insert remove 
    #max order 2 
    #add extend

        