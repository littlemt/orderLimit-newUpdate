#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:19:13 2022

@author: littlemt
"""

import configparser
import numpy as np

rng=np.random.default_rng()



a=np.linspace(1, 10,10)
a=np.reshape(a, (2,5))


print(np.where(2==a)[0])
'''


j=rng.integers(0,12)
print(j,'order')
a=np.linspace(0, j,j+1)
c=np.zeros(14)
c[:j+1]=a
c[j+1]=j+1
e=c


l=rng.integers(0,j+1)
for i in range(100):
    
    
    
    print(l,'random')
    
    
    dum=np.random.uniform(c[l],c[l+1])
    dum2=np.random.uniform(dum,c[l+1])
                          
                          
    b=np.array([c[l],dum,dum2,c[l+1]])
    d=np.array([c[l],c[l+1]])
    
    
    print(c,'ins')
    length=len(b)
    index1=l
    index2=j+2
    #print(c[index1+length:index2+length],c[index1:index2])
    c[index1+length-2:index2+length-2]=c[index1:index2]
    c[index1+1:index1+length-1]=b[1:-1]
    i+=1
    print(c,'rem')
    
    
    
    length=len(d)
    c[index1:index1+length]=d
    c[index1+length:index2]=c[index1+length+2:index2+2]
    c[index2:]=0
    
    print(c)
    '''