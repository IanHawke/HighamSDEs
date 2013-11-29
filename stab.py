# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 22:49:20 2013

stab

This is a direct re-implementation of Des Higham's SDE scripts. 
Mean-square and asymptotic stability test for Euler-Maruyama

SDE is
    dX = \lambda X dt + \mu X dW,   X(0) = X_0
with \lambda = -3, \mu = \sqrt{3}, X_0 = 1.

@author: ih3
"""

import numpy as np
import matplotlib.pyplot as plt

T = 20.0
M = 50000
X0 = 1.0
lam = -3.0
mu = np.sqrt(3.0)

dt = 0.25
N = 80
plt.subplot(2,1,1)
for k in range(3):
    Dt = 2**k*dt
    L = N/2**k
    Xms = np.zeros(L)
    Xtemp = X0 * np.ones(M)
    t = np.linspace(0,T,L)
    for j in range(L):
        Winc = np.sqrt(Dt)*np.random.randn(M)
        Xtemp += (Dt*lam + mu*Winc)*Xtemp
        Xms[j] = np.mean(Xtemp**2)
    plt.semilogy(t,Xms)
plt.legend(('\Delta t = 1/4', '\Delta t = 1/2', '\Delta t = 1'))
plt.title('Mean-Square: \lambda = -3, \mu = \sqrt{3}')
plt.ylabel('E[X^2]')

plt.subplot(2,1,2)
T = 500.0
lam = 0.5
mu = np.sqrt(6)
N = 2000
for k in range(3):
    Dt = 2**k*dt
    L = N/2**k
    Xemabs = np.zeros(L)
    Xtemp = X0
    t = np.linspace(0,T,L)
    for j in range(L):
        Winc = np.sqrt(Dt)*np.random.randn(1)
        Xtemp += (Dt*lam + mu*Winc)*Xtemp
        Xemabs[j] = np.abs(Xtemp)
    plt.semilogy(t,Xemabs)
plt.legend(('\Delta t = 1/4', '\Delta t = 1/2', '\Delta t = 1'))
plt.title('Single Path: \lambda = 1/2, \mu = \sqrt{6}')
plt.ylabel('|X|')


plt.show()
