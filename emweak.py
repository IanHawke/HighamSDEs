# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 20:11:14 2013

emweak

This is a direct re-implementation of Des Higham's SDE scripts. 
Testing weak convergence of the Euler-Maruyama method applied to a linear SDE

SDE is
    dX = \lambda X dt + \mu X dW,   X(0) = X_0
with \lambda = 2, \mu = 1, X_0 = 1

Uses 2^k steps where k = 5, ..., 9.

Weak convergence is | E(X_L) - E(X(T)) |

@author: ih3
"""

import numpy as np

EM_weak = False

lam = 2.0
mu = 0.1
X0 =1.0

T = 1.0
N = 2**9
dt = T / N

M = 50000 # Number of paths sampled at each resolution

Xem = np.zeros(5)
Dtvals = np.zeros(5)

for p in range(5):
    R = 2**p
    Dt = R*dt
    Dtvals[p] = Dt
    L = N / R
    Xtemp = X0 * np.ones(M)
    for j in range(L):
        if EM_weak:
            Winc = np.sqrt(Dt)*np.sign(np.random.randn(M)) # Weak E-M algorithm
        else:
            Winc = np.sqrt(Dt)*np.random.randn(M)
        Xtemp += (Dt * lam + mu*Winc)*Xtemp
    Xem[p] = np.mean(Xtemp)
Xerr = np.abs(Xem - np.exp(lam))

# Measure the convergence rate
p = np.polyfit(np.log(Dtvals),np.log(Xerr),1)
print "Measured convergence rate is ", p[0]

import matplotlib.pyplot as plt
plt.loglog(Dtvals,Xerr,'b*-')
plt.loglog(Dtvals,np.exp(p[1])*Dtvals,'r--')
plt.xlabel('\Delta t')
plt.ylabel('|E(X(t)) - Sample average of X_L|')
plt.show()

