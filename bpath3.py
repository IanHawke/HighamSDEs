# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:37:59 2013

bpath3

This is a direct re-implementation of Des Higham's SDE scripts. 
Function along a Brownian path

@author: ih3
"""

import numpy as np

T = 1.0 # End time
N = 500 # Number of steps
dt = T / N

M = 1000 # Will do M paths simultaneously

dW = np.sqrt(dt)*np.random.randn(M,N)
W = np.cumsum(dW,1)
t = np.linspace(0,T,N)
U = np.exp(t + 0.5*W)
Umean = np.mean(U,0)

print "Average error is ", np.linalg.norm((Umean - np.exp(9*t/8)),np.inf)

import matplotlib.pyplot as plt

plt.plot(t,Umean,'b-')
for i in range(5):
    plt.plot(t,U[i,:],'r--')
plt.xlabel('t')
plt.ylabel('U(t)')
plt.legend(('Mean of all paths','Individual paths'))
plt.show()
