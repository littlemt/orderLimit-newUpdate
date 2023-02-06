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
import os
import shutil as st

from matplotlib import rc 
rc('text', usetex=True)





#need to define a function that takes makes a new random generator 
def plot0(hist,p,mu,directory='./',m=1):
    
    
    x=hist[:,0]
    y=hist[:,1]
    yerr=hist[:,2]
    
             
    mpl.xlabel(r'$\tau$')
    mpl.ylabel('log[-G(p=0,tau)]')
    mpl.title(r'$\mu$='+str(mu))
    mpl.errorbar(x,np.log(-hist[:,1]),yerr=yerr/y ,fmt='o',label='Data')
    mpl.plot(x,-x*(p/2/m-mu),color='orange',zorder=3,label='Exact')
    m,b=np.polyfit(x,np.log(-y),deg=1)
    mpl.plot(x,m*x+b,color='red',zorder=2,label='regression')
    mpl.title('reg line: '+'y='+str(round(m,5))+'x+'+str(round(b,5)))
    mpl.savefig(directory+'tauvsLogG0_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    mpl.legend()
    
    mpl.show()
    
    mpl.errorbar(x,-hist[:,1] ,yerr=hist[:,2] ,fmt='o',label='Data')
    mpl.plot(hist[:,0],np.exp(hist[:,0]*mu),color='red',zorder=2)
    mpl.xlabel(r'$\tau$')
    mpl.ylim(0,1)
    mpl.ylabel('G')
    mpl.savefig(directory+'tauvsG0_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    mpl.show()
    
    mpl.errorbar(x,-hist[:,1] -np.exp(-x*(p/2/m-mu)),yerr=hist[:,2] ,fmt='o',label='Data')
    mpl.xlabel(r'$\tau$')
    mpl.ylabel('G')
    mpl.savefig(directory+'tauvsG0_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    mpl.show()
    


def plot1(hist,count,order,p,mu,directory='./',m=1):
    
    
    x=hist[:,0]
    y=hist[:,1]
    yerr=hist[:,2]
             
    mpl.xlabel('tau')
    mpl.ylabel('log[-G(p=0,tau)]')
    mpl.title('mu='+str(mu))
    mpl.errorbar(x,np.log(-hist[:,1]),yerr=yerr/y ,fmt='o',label='Data')
    mpl.plot(x,np.log(np.exp(-(p**2/(2*m)-mu)*x)-fpc.firstOrderSolution(x, mu)),color='orange',zorder=2,label='Exact')
    m,b=np.polyfit(x[int(.25*len(x)):],np.log(-y[int(.25*len(x)):]),deg=1)
    mpl.plot(x[int(.25*len(x)):],m*x[int(.25*len(x)):]+b,color='red',zorder=3,label='regression')
    mpl.title('reg line: '+'y='+str(round(m,5))+'x+'+str(round(b,5)))
    mpl.legend()
    mpl.savefig(directory+'tauvsLogG1_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    
    mpl.show()
    
    mpl.xlabel('order')
    mpl.ylabel('% of time')
    mpl.bar(np.arange(len(order)),order/sum(order))
    mpl.title('delta MC order 0 ='+str(count[0]))
    mpl.savefig(directory+'ordPlot'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    mpl.show()
    
    mpl.bar(['insert %','remove %'],[count[3]*100/count[4],count[5]*100/count[6]])
    mpl.savefig(directory+'accProb_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    mpl.show()#fix the monte carlo time in this
    
    mpl.errorbar(x,-hist[:,1] -np.exp(-x*(p/2/m-mu))-fpc.firstOrderSolution(x,mu),yerr=hist[:,2] ,fmt='o',label='Data')
    mpl.xlabel(r'$\tau$')
    mpl.ylabel('G')
    mpl.savefig(directory+'tauvsG0_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    mpl.show()
    
    
def loop(seed):
    #this code take in the random seed then parse the param file then run the loop
    
    
    
    config=configparser.ConfigParser()
    
    config.read('param.ini')
    
    tauMax=float(config.get('section_a','tauMax'))
    runTime=int(config.get('section_a','runTime'))
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
        
        mass=float(config.get('section_a','mass'))
        omega=float(config.get('section_a','omega'))
        pExt=float(config.get('section_a','exMomentum'))
        mu=float(config.get('section_a','mu'))
        maxOrder=int(config.get('section_a','maxOrder'))
        
        histArray=np.zeros((bins,3))
        histArray[:,0]=result[0][0][:,0]
        noZero=np.zeros(noThread)
        count=np.zeros((12,noThread))
        order=np.zeros((maxOrder+1,noThread))
        countAvg=np.zeros(12)
        orderAvg=np.zeros(maxOrder+1)
        
        
       
        histR=np.ndarray((bins,noThread))
        for i in range(noThread):
            noZero[i]=result[i][1]
            histR[:,i]=fpc.calc(result[i][0][:,1],histArray[-1,0],histArray[1,0],pExt,mu,result[i][1])
            count[:,i]=result[i][2]
            order[:,i]=result[i][3]
            
    
        for i in range(maxOrder+1):
            orderAvg[i]=np.average(order[i])
            
        for i in range(12):
            countAvg[i]=np.average(count[i])
        
        for i in range(bins):
            histArray[i,1]=np.average(histR[i])
            histArray[i,2]=np.std(histR[i])/noThread**.5 
            
        directory='./Plots/O'+config.get('section_a','maxOrder')+'_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_T'+config.get('section_a','tauMax')+'/'
        os.mkdir(directory)
        
        np.savetxt(directory+'avgData',histArray,delimiter=',')
        np.savetxt(directory+'binData',histR,delimiter=',')
        np.savetxt(directory+'zeros', noZero,delimiter=',')
        np.savetxt(directory+'counts',count,delimiter=',')
        
        if int(config.get('section_a','maxOrder'))==0:
            plot0(histArray,pExt,mu,directory=directory)
        elif int(config.get('section_a','maxOrder'))==1:
            plot1(histArray,countAvg,orderAvg,pExt,mu,directory=directory)
        else:
            
            mpl.xlabel(r'$\tau$', fontsize=18)
            mpl.ylabel(r'$\log[-G(p=0,\tau)]$')
            mpl.title(r'$\mu=$'+str(mu))
            mpl.errorbar(histArray[:,0],np.log(-histArray[:,1]),yerr=histArray[:,2] ,fmt='o')
            mpl.savefig(directory+'tauvsLogG_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
            
            mpl.show()
            
            mpl.xlabel('order')
            mpl.ylabel('% of time')
            mpl.bar(np.arange(len(order)),order/sum(order))
            mpl.title('delta MC order 0 ='+str(count[0]))
            mpl.savefig(directory+'orderHist'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
            mpl.show()
            
            mpl.bar(['insert %','remove %','swap %'],[count[3]*100/count[4],count[5]*100/count[6],count[7]*100/count[8]])
            mpl.savefig(directory+'accProb'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+str(pExt)+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
            mpl.show()
        
        
        
            
        
        
        
        
        #run the function deffined above per thread combine to large data set
        #take outup then run the calculation step in paralell  
        



