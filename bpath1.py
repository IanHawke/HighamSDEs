# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:37:59 2013

bpath1

This is a direct re-implementation of Des Higham's SDE scripts. First up
is the Brownian path simulation

@author: ih3
"""

import numpy as np

T = 1.0 # End time
N = 500 # Number of steps
dt = T / N

W = np.zeros(N)
dW = np.zeros(N)
t = np.linspace(0,T,N)

dW[0] = np.random.randn()
W[0] = dW[0]

for j in range(1,N):
    dW[j] = np.sqrt(dt) * np.random.randn()
    W[j] = W[j-1] + dW[j]
    
import matplotlib.pyplot as plt

plt.plot(t,W)
plt.xlabel('t')
plt.ylabel('W(t)')
plt.show()
