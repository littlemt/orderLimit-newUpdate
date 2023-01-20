#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:11:36 2022

@author: littlemt
"""

from multiprocessing import Pool
import fpc_orderLimit as fpc
import configparser
import numpy as np
from matplotlib import pyplot as mpl





#need to define a function that takes makes a new random generator 

def loop(seed):
    #this code take in the random seed then parse the param file then run the loop
    
    
    
    config=configparser.ConfigParser()
    
    config.read('param.ini')
    
    tauMax=float(config.get('section_a','tauMax'))
    runTime=float(config.get('section_a','runTime'))
    upProbs=np.int16(config.get('section_a','updateProb').split(','))
    pExt=float(config.get('section_a','exMomentum'))
    mu=float(config.get('section_a','mu'))
    alpha=float(config.get('section_a','alpha'))
    mass=float(config.get('section_a','mass'))
    omega=float(config.get('section_a','omega'))
    orderMax=int(config.get('section_a','maxOrder'))
    thermal=int(config.get('section_a','thermal'))
    step=int(config.get('section_a','step'))
    bins=int(config.get('section_a','bins'))
    
    #need to define all paramaters from file
    
    #run the code for the loop to collect the data
    
    data=fpc.first_order(tauMax,runTime,upProbs,pExt,mu,alpha,orderMax,thermal,step,seed,bins=bins,omega=omega,m=mass,debug=0)
    
    
    #fcp.saveData(data,'/home/littlemt/Documents/plots_data/',tauMax,runTime,upProbs,pExt,mu,alpha,orderMax,thermal,step,bins,seed)
    return data
    #maybe just return the data



if __name__ =='__main__':
    
    config=configparser.ConfigParser()

    config.read('param.ini')
    seed=int(config.get('section_b','seed'))
    noThread=int(config.get('section_b','noThread'))
    
    
    rng=np.random.default_rng(seed)
    
    rand=rng.integers(0,int(1E10),noThread)
    #gen random numbers same dim as number of threds
    #run the parallel using the random num as seeds
    bins=int(config.get('section_a','bins'))
    
    with Pool(processes=(noThread)) as p:
        
        result=p.map(loop,rand)
        
        
        histArray=np.zeros((bins,3))
        histArray[:,0]=result[0][0][:,0]
        noZero=0
        count=np.zeros(12)
        order=np.zeros(bins+1)
        for i in result:
            histArray[:,1]+=i[0][:,1]
            noZero+=i[1]
            count+=i[2]
            order+=i[3]
            
            
        
        mass=float(config.get('section_a','mass'))
        omega=float(config.get('section_a','omega'))
        pExt=float(config.get('section_a','exMomentum'))
        mu=float(config.get('section_a','mu'))
    
        mpl.xlabel('tau')
        mpl.ylabel('log[-G(p=0,tau)]')
        mpl.title('mu='+str(mu))
        mpl.scatter(histArray[:,0],np.log(-fpc.calc(histArray,pExt,mu,noZero)))
        
        
        mpl.show()
        
        mpl.xlabel('order')
        mpl.ylabel('% of time')
        mpl.bar(np.arange(len(order)),order/sum(order))
        mpl.title('delta MC order 0 ='+str(count[0]))
        
        mpl.show()
        
        mpl.bar(['insert %','remove %','swap %'],[count[3]*100/count[4],count[5]*100/count[6],count[7]*100]/count[8])
        mpl.show()
        
        
        
            
        
        
        
        
        #run the function deffined above per thread combine to large data set
        #take outup then run the calculation step in paralell  
        



