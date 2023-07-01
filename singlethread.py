#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 09:09:46 2023

@author: littlemt
"""

import fpc_orderLimit as fpc
import plotting as fplot
import configparser
import numpy as np
import os






#need to define a function that takes makes a new random generator 

    
def loop(seed):
    #this code take in the random seed then parse the param file then run the loop
    
    
    
    config=configparser.ConfigParser()
    
    config.read('param.ini')
    
    tauMax=float(config.get('section_a','tauMax'))
    runTime=int(config.get('section_a','runTime'))
    upProbs=np.float64(config.get('section_a','updateProb').split(','))
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

def jackknife(a):
    '''
    this function runs jackknife satistics on input array a

    Parameters
    ----------
    a : array
        input array.

    Returns
    -------
    None.

    '''
    dim=np.shape(a)
    dummyArray=np.ndarray(dim)

    for i in range(dim[0]):
        dummyArray[i]=np.average(jkSlice(i,a),axis=0)
        
        
    return dummyArray

def jkSlice(index,array):
    #could probobly do this destructivly dont know if this is more efficient 
    #or not but it should not matter 
    dim=list(np.shape(array))
    dim[0]=dim[0]-1
    dummyArray=np.zeros(tuple(dim))
    
    dummyArray[:index]=array[:index]
    dummyArray[index:]=array[index+1:]
    
    return dummyArray
    
    

if __name__ =='__main__':
    
    config=configparser.ConfigParser()

    config.read('param.ini')
    seed=int(config.get('section_b','seed'))
    noHist=int(config.get('section_b','noHist'))
    
    rng=np.random.default_rng(seed)
    
    rand=rng.integers(0,int(1E10))
    #gen random numbers same dim as number of threds
    #run the parallel using the random num as seeds
    bins=int(config.get('section_a','bins'))
    

    mass=float(config.get('section_a','mass'))
    omega=float(config.get('section_a','omega'))
    pExt=float(config.get('section_a','exMomentum'))
    mu=float(config.get('section_a','mu'))
    maxOrder=int(config.get('section_a','maxOrder'))
    alpha=float(config.get('section_a','alpha'))
    
    
    result=np.ndarray((noHist,4),dtype=object)
    
    for i in range(noHist):
        result[i]=loop(rng.integers(0,int(1E10)))
     
    histArray=np.zeros((bins,3))
    histArray[:,0]=result[0][0][:,0]
    noZero=np.zeros(noHist)
    count=np.zeros((noHist,13))
    order=np.zeros((noHist,maxOrder+1))
    countAvg=np.zeros(12)
    orderAvg=np.zeros(maxOrder+1)
    
    
   
    histR=np.ndarray((noHist,bins))
    
    for i in range(noHist):
        #this loop unpacks the output from the pool
        noZero[i]=result[i][1]
        histR[i,:]=result[i][0][:,1]
        count[i,:]=result[i][2]
        order[i,:]=result[i][3]
        
    jkArray=jackknife(histR)
    jkZero=jackknife(noZero)
    jkCount=jackknife(count)
    
    jkHist=np.ndarray(np.shape(jkArray))
    for i in range(noHist):
        jkHist[i]=fpc.calc(jkArray[i],histArray[-1,0],histArray[0,0]*2,pExt,mu,jkZero[i])
        

    #,histArray[-1,0],histArray[0,0]*2,pExt,mu,result[i][1])

    
    orderAvg=np.average(order,axis=0)
    countAvg=np.average(count,axis=0)
    histArray[:,1]=np.average(jkHist,axis=0)
    
    #is this error estimation proper or do you have to do it differnetly with jackknife?
    #there was some thing about it in the notes by peter young but dont know what they mean?
    #not excatly sure what eq 53 means?
    
    histArray[:,2]=np.std(jkHist,axis=0)*(noHist-1)**.5/noHist**.5 
        

    directory='./Plots/Single_noRe_O'+config.get('section_a','maxOrder')+'_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_T'+config.get('section_a','tauMax')+'/'

    
    try:
        os.mkdir(directory)
    except:
        print('directory already there')
        
    
    np.savetxt(directory+'avgData_Seed'+config.get('section_b','seed'),histArray,delimiter=',')
    np.savetxt(directory+'binData_Seed'+config.get('section_b','seed'),histR,delimiter=',')
    np.savetxt(directory+'zeros_Seed'+config.get('section_b','seed'), noZero,delimiter=',')
    np.savetxt(directory+'counts_Seed'+config.get('section_b','seed'),count,delimiter=',')
    
    if int(config.get('section_b','plot'))==1:
        if int(config.get('section_a','maxOrder'))==0:
            fplot.plot0(histArray,pExt,mu,directory=directory)
        elif int(config.get('section_a','maxOrder'))==1:
            fplot.plot1(histArray,countAvg,orderAvg,pExt,mu,alpha,directory=directory)
        else:
            fplot.plot(histArray,countAvg,orderAvg,pExt,mu,directory=directory)
        
    
    
        
    
    
    
    
    #run the function deffined above per thread combine to large data set
    #take outup then run the calculation step in paralell  
    
