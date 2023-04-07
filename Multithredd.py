#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:11:36 2022

@author: littlemt
"""

from multiprocessing import Pool
import fpc_orderLimit as fpc
import plotting as fplot
import configparser
import numpy as np
import os


from matplotlib import rc 
rc('text', usetex=True)





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
        
        mass=float(config.get('section_a','mass'))
        omega=float(config.get('section_a','omega'))
        pExt=float(config.get('section_a','exMomentum'))
        mu=float(config.get('section_a','mu'))
        maxOrder=int(config.get('section_a','maxOrder'))
        alpha=float(config.get('section_a','alpha'))
         
        histArray=np.zeros((bins,3))
        histArray[:,0]=result[0][0][:,0]
        noZero=np.zeros(noThread)
        count=np.zeros((13,noThread))
        order=np.zeros((maxOrder+1,noThread))
        countAvg=np.zeros(12)
        orderAvg=np.zeros(maxOrder+1)
        
        
       
        histR=np.ndarray((bins,noThread))
        
        for i in range(noThread):
            #this loop unpacks the output from the pool
            noZero[i]=result[i][1]
            histR[:,i]=fpc.calc(result[i][0][:,1],histArray[-1,0],histArray[0,0]*2,pExt,mu,result[i][1])
            count[:,i]=result[i][2]
            order[:,i]=result[i][3]
            
    
        
        orderAvg=np.average(order,axis=1)
        countAvg=np.average(count,axis=1)
        histArray[:,1]=np.average(histR,axis=1)
        histArray[:,2]=np.std(histR,axis=1)/noThread**.5 
            

        directory='./Plots/noRe_O'+config.get('section_a','maxOrder')+'_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_T'+config.get('section_a','tauMax')+'/'

        
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
        



