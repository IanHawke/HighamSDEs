# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 19:14:36 2013

stint

This is a direct re-implementation of Des Higham's SDE scripts. 
Stochastic integrals, comparing Ito and Stratonovich cases.
Results are different values for \int W dW


@author: ih3
"""

import numpy as np

T = 1.0
N = 500
dt = T / N

dW = np.sqrt(dt) * np.random.randn(N)
W = np.cumsum(dW)

tmpW = np.zeros(N)
tmpW[1:] = W[:-1]

ito = sum(tmpW*dW)
strat = sum((0.5*(tmpW+W) + 0.5*np.sqrt(dt)*np.random.randn(N))*dW)

print "The Ito integral is ", ito, " with error ", np.abs(ito - 0.5*(W[-1]**2-T))
print "The Stratonovich integral is ", strat, " with error ", np.abs(strat - 0.5*W[-1]**2)
