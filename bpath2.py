# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:37:59 2013

bpath2

This is a direct re-implementation of Des Higham's SDE scripts. First up
is the Brownian path simulation

@author: ih3
"""

import numpy as np

T = 1.0 # End time
N = 500 # Number of steps
dt = T / N

dW = np.sqrt(dt)*np.random.randn(N)
W = np.cumsum(dW)
t = np.linspace(0,T,N)

    
import matplotlib.pyplot as plt

plt.plot(t,W)
plt.xlabel('t')
plt.ylabel('W(t)')
plt.show()
