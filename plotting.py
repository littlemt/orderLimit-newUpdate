#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 20:28:31 2023

@author: littlemt
"""
import matplotlib.pyplot as mpl
from matplotlib import rc 
rc('text', usetex=True)
import numpy as np
import fpc_orderLimit as fpc
import configparser
config=configparser.ConfigParser()

from scipy.special import erf

def firstOrderSolution(tau,mu,alpha,omega=1,m=1):
    return -2*alpha*2**.5*np.pi*np.exp(tau*(mu-omega))*m**.5*(2*(omega*tau)**.5+np.exp(omega*tau)*np.pi**.5*(2*omega*tau-1)*erf((omega*tau)**.5))/(32**.5*omega**(1.5)*np.pi**1.5)




def plot1(hist,count,order,p,mu,alpha,directory='./',m=1):
    config.read('param.ini')
    
    x=hist[:,0]-hist[0,0]
    y=hist[:,1]
    yerr=hist[:,2]
             
    mpl.xlabel(r'$\tau$')
    mpl.ylabel('log[-G(p=0,tau)]')
    mpl.title(r'$mu=$'+str(mu))
    mpl.errorbar(x,np.log(-hist[:,1]),yerr=np.abs(yerr/y) ,fmt='o',label='Data')
    mpl.plot(x,np.log(np.exp(-(p**2/(2*m)-mu)*x)-firstOrderSolution(x, mu,alpha)),color='orange',zorder=2,label='Exact')

    mpl.legend()
    mpl.xlim(x[0],x[-1])
    mpl.savefig(directory+'tauvsLogG1_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    
    mpl.show()
    
    mpl.xlabel('order')
    mpl.ylabel('% of time')
    mpl.bar(np.arange(len(order)),order/sum(order))
    mpl.title('delta MC order 0 ='+str(count[0]))
    mpl.savefig(directory+'ordPlot'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    mpl.show()

    
    mpl.errorbar(x,np.log(-hist[:,1]) -np.log(np.exp(-x*(p/2/m-mu))+firstOrderSolution(x,mu,alpha)),yerr=hist[:,2] ,fmt='o',label='Data')
    mpl.xlabel(r'$\tau$')
    mpl.ylabel('G')
    mpl.xlim(x[0],x[-1])
    mpl.plot(x,0*x,zorder=2)
    mpl.savefig(directory+'tauvsG-acc0_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    mpl.show()
    
    mpl.errorbar(x,-hist[:,1] -np.exp(-x*(p/2/m-mu))+firstOrderSolution(x,mu,alpha),yerr=hist[:,2] ,fmt='o',label='Data')
    mpl.xlabel(r'$\tau$')
    mpl.ylabel('G')
    mpl.xlim(x[0],x[-1])
    mpl.plot(x,0*x,zorder=2)
    mpl.savefig(directory+'tauvsG-acc0_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    mpl.show()
    
    print(np.average(-hist[:,1] -np.exp(-x*(p/2/m-mu))+firstOrderSolution(x,mu,alpha)),np.std(-hist[:,1] -np.exp(-x*(p/2/m-mu))-fpc.firstOrderSolution(x,mu,alpha)))
    
def plot0(hist,p,mu,directory='./',m=1):
    config.read('param.ini')
    
    x=hist[:,0]
    y=hist[:,1]
    yerr=hist[:,2]
    
             
    mpl.xlabel(r'$\tau$')
    mpl.ylabel(r'$\log[-G(p=0,\tau)]$')
    mpl.title(r'$\mu$='+str(mu))
    mpl.errorbar(x,np.log(-hist[:,1]),yerr=yerr/abs(y) ,fmt='o',label='Data')
    mpl.plot(x,-(p**2/2/m-mu)*x,color='orange',zorder=3,label='Exact')
    mpl.xlim(x[0],x[-1])
    mpl.legend()
    mpl.savefig(directory+'tauvsLogG0_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    
    mpl.show()
    
    mpl.errorbar(x,-hist[:,1] ,yerr=hist[:,2] ,fmt='o',label='Data')
    mpl.plot(hist[:,0],(np.exp(-x*(p/2/m-mu))),color='red',zorder=2)
    mpl.xlabel(r'$\tau$')
    mpl.ylim(0,1)
    mpl.ylabel('G')
    mpl.xlim(x[0],x[-1])
    mpl.savefig(directory+'tauvsG0_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    mpl.show()
    
    mpl.errorbar(x,-hist[:,1] -np.exp(-x*(p/2/m-mu)),yerr=hist[:,2] ,fmt='o',label='Data',zorder=1)
    mpl.xlabel(r'$\tau$')
    mpl.plot(x,0*x,zorder=2)
    mpl.xlim(x[0],x[-1])
    mpl.savefig(directory+'tauvsG0-acc_m'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    mpl.show()
    
    
    
def plot(hist,count,order,p,mu,directory='./',m=1):
    config.read('param.ini')
    x=hist[:,0]
    y=hist[:,1]
    yerr=hist[:,2]
    
   
    mpl.xlabel(r'$\tau$', fontsize=18)
    mpl.ylabel(r'$\log[-G(p=0,\tau)]$')
    mpl.title(r'$\mu=$'+str(mu))
    mpl.errorbar(x,np.log(-y),yerr=yerr ,fmt='o')
    mpl.xlim(x[0],x[-1])
    mpl.ylim(-20,0)
    mpl.savefig(directory+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    mpl.show()

    
    mpl.xlabel('order')
    mpl.ylabel('% of time')
    mpl.bar(np.arange(len(order)),order/sum(order))
    mpl.title('delta MC order 0 ='+str(count[0]))
    mpl.savefig(directory+'orderHist'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    mpl.show()
    
    # mpl.bar(['insert %','remove %','swap %'],[count[3]*100/count[4],count[5]*100/count[6],count[7]*100/count[8]])
    # mpl.savefig(directory+'accProb'+str(mu)+'_P='+config.get('section_a','updateProb')+'_p'+config.get('section_a','exMomentum')+'_a'+config.get('section_a','alpha')+'_rt'+config.get('section_a','runTime')+'_O'+config.get('section_a','maxOrder')+'.pdf' )
    
    # mpl.show()
    