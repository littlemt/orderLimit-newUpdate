#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:11:36 2022

@author: littlemt
"""

from multiprocessing import Pool
import fpc_orderLimit as fcp
import configparser
import numpy as np





#need to define a function that takes makes a new random generator 

def loop(n):
    #this code take in the random seed then parse the param file then run the loop
    
    rng=np.random.default_rng(n)
    config=configparser.ConfigParser()

    config.read('param.ini')
    
    tauMax=float(config.get('section_a','tauMax'))
    runTime=float(config.get('section_a','runTime'))
    upProbs=np.int16(config.get('section_a','updateProb').split(','))
    pExt=float(config.get('section_a','exMomentum'))
    mu=float(config.get('section_a','mu'))
    k=float(config.get('section_a','k'))
    alpha=float(config.get('section_a','alpha'))
    mass=float(config.get('section_a','mass'))
    omega=float(config.get('section_a','omega'))
    orderMax=int(config.get('section_a','maxOrder'))
    
    
    
    #need to define all paramaters from file
    
    #run the code for the loop to collect the data
    
    data=fcp.full(tauMax,runTime,upProbs,pExt,mu,k,alpha,orderMax,mass)
    
    return data
    #maybe just return the data
    

if __name__ =='__main__':
    
    config=configparser.ConfigParser()

    config.read('param.ini')
    seed=float(config.get('section_b','seed'))
    noThread=float(config.get('section_b','noThread'))
    
    
    rng=np.random.default_rng(seed)
    
    rng.uniform(0,1000,noThread)
    #gen random numbers same dim as number of threds
    #run the parallel using the random num as seeds
    
    
    with Pool as p:
        
        p.map(loop,)
        
        #run the function deffined above per thread combine to large data set
        #take outup then run the calculation step in paralell  
        



