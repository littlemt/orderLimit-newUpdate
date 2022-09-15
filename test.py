#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:19:13 2022

@author: littlemt
"""

import configparser
import numpy as np

a=np.linspace(1, 10,10)
b=np.reshape(a, (5,2))
c=np.where(b==10)

for i in range(2,4):
    print(i)