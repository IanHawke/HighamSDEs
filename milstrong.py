# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 20:11:14 2013

milstrong

This is a direct re-implementation of Des Higham's SDE scripts. 
Testing strong convergence of the Milstein method applied to the SDE
    dX = r X (K - X) dt + \beta X dW,   X(0) = X_0
with r = 2, K = 1, \beta = 0.25, X_0 = 0.5

Uses 2^k steps where k = 4, ..., 7 against a reference solution with 2^11 steps.

Strong convergence is E | X_L - X(T) |

@author: ih3
"""

import numpy as np


r = 2.0
K = 1.0
beta = 0.25
X0 =0.5

T = 1.0
N = 2**11
dt = T / N

M = 500 # Number of paths sampled at each resolution

dW = np.sqrt(dt)*np.random.randn(M,N)
Xmil = np.zeros((M,5))

Dtvals = np.zeros(5)

R = [1, 16, 32, 64, 128]
for p in range(5):
    Dt = R[p]*dt
    Dtvals[p] = Dt
    L = N / R[p]
    Xtemp = X0 * np.ones(M)
    for j in range(L):
        Winc = np.sum(dW[:,R[p]*j:R[p]*(j+1)],1)
        Xtemp += Dt*r*Xtemp*(K-Xtemp) + beta*Xtemp*Winc + 0.5*beta**2*Xtemp*(Winc**2 - Dt)
    Xmil[:,p] = Xtemp

Xref = Xmil[:,0]
Xerr = np.zeros(4)
for i in range(1,5):
    Xerr[i-1] = np.mean(np.abs(Xmil[:,i] - Xref))

# Measure the convergence rate
p = np.polyfit(np.log(Dtvals[1:]),np.log(Xerr),1)
print "Measured convergence rate is ", p[0]

import matplotlib.pyplot as plt
plt.loglog(Dtvals[1:],Xerr,'b*-')
plt.loglog(Dtvals[1:],np.exp(p[1])*Dtvals[1:],'r--')
plt.xlabel('\Delta t')
plt.ylabel('Sample average of |X(t) - X_L|')
plt.show()

