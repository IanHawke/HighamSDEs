# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 20:11:14 2013

emstrong

This is a direct re-implementation of Des Higham's SDE scripts. 
Testing strong convergence of the Euler-Maruyama method applied to a linear SDE

SDE is
    dX = \lambda X dt + \mu X dW,   X(0) = X_0
with \lambda = 2, \mu = 1, X_0 = 1

Uses 2^k steps where k = 5, ..., 9.

Strong convergence is E | X_L - X(T) |

@author: ih3
"""

import numpy as np

lam = 2.0
mu = 1.0
X0 =1.0

T = 1.0
N = 2**9
dt = T / N

M = 1000 # Number of paths sampled at each resolution

Xerr = np.zeros((M,5))

for s in range(M):
    dW = np.sqrt(dt) * np.random.randn(N)
    W = np.cumsum(dW)
    Xtrue = X0 * np.exp((lam - 0.5*mu**2)*T + mu*W[-1])
    for p in range(5):
        R = 2**p
        Dt = R*dt
        L = N / R
        Xtemp = X0
        for j in range(L):
            Winc = sum(dW[R*j:R*(j+1)])
            Xtemp += (Dt*lam + mu*Winc)*Xtemp
        Xerr[s,p] = np.abs(Xtemp - Xtrue)

Dtvals = Dt * 2**(np.array(range(5)))

# Measure the convergence rate
p = np.polyfit(np.log(Dtvals),np.log(np.mean(Xerr,0)),1)
print "Measured convergence rate is ", p[0]

import matplotlib.pyplot as plt
plt.loglog(Dtvals,np.mean(Xerr,0),'b*-')
plt.loglog(Dtvals,np.exp(p[1])*np.sqrt(Dtvals),'r--')
plt.xlabel('\Delta t')
plt.ylabel('Sample average of |X(t) - X_L|')
plt.show()

