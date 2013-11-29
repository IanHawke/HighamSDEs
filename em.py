# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 19:14:36 2013

em

This is a direct re-implementation of Des Higham's SDE scripts. 
The Euler-Maruyama method applied to a linear SDE

SDE is
    dX = \lambda X dt + \mu X dW,   X(0) = X_0
with \lambda = 2, \mu = 1, X_0 = 1


@author: ih3
"""

import numpy as np

lam = 2.0
mu = 1.0
X0 =1.0

T = 1.0
N = 2**8
dt = T / N

dW = np.sqrt(dt) * np.random.randn(N)
W = np.cumsum(dW)
t = np.linspace(0,T,N)

Xtrue = X0 * np.exp((lam - 0.5*mu**2)*t + mu*W)

R = 4
Dt = R * dt
L = N / R
Xem = np.zeros(L)
Xem[0] = X0

Xtemp = X0
for j in range(L):
    Winc = sum(dW[R*j:R*(j+1)])
    Xtemp += (Dt*lam + mu*Winc)*Xtemp
    Xem[j] = Xtemp

import matplotlib.pyplot as plt

tt = np.linspace(0,T,L)
plt.plot(t,Xtrue,'m-')
plt.plot(tt,Xem,'r--*')
plt.xlabel('t')
plt.ylabel('X')
plt.show()

print "Error at the endpoint is ", np.abs(Xem[-1] - Xtrue[-1])