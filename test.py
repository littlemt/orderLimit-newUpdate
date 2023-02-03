#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:19:13 2022

@author: littlemt
"""

import configparser
import numpy as np
import matplotlib.pyplot as plt



rng=np.random.default_rng()

x=np.linspace(1,20,20)
y=x**1.5
err=rng.uniform(0,2,20)**2

plt.errorbar(x,y,yerr=err,fmt='o')
plt.savefig('test.pdf')