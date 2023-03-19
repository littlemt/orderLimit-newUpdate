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
from scipy.integrate import nquad 

def func(theta,k1,k2,tau4,tau3,tau2,tau1):
    return 1