#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:19:13 2022

@author: littlemt
"""

import configparser
import numpy as np
import matplotlib.pyplot as plt
from  math import floor 



rng=np.random.default_rng()

tauMax=15
bins=29
deltaTau=tauMax/(bins+1)
print(deltaTau)

x=np.linspace(deltaTau,tauMax-deltaTau,bins)
print(x)



hist=np.zeros((bins,2))
hist[:,0]=x
h=[]
for t in np.arange(0,15,.25):
    
    h.append(t)
    hist[floor(t/deltaTau),1]+=1
    print(t,int(t/deltaTau))
    


plt.scatter(x,hist[:,1])
