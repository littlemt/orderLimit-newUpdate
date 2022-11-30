#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:19:13 2022

@author: littlemt
"""

import configparser
import numpy as np



rng=np.random.default_rng()

def uni():
    return rng.uniform()



def norm(array):
    return np.linalg.norm(array)
normVec=np.vectorize(norm)

def normVec1(array):
    dumArray=np.zeros(len(array))
    for i in range(len(array)):
        dumArray[i]=np.linalg.norm(array[i])
    return dumArray

a=np.array([[1,2,3],[1,2,3],[1,2,3]])
b=np.array([[1,2,3],[1,2,3],[1,2,3]],dtype=object)

print(norm([1,2,3]))
print(normVec(b))
print(normVec1(a))