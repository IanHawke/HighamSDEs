import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

def mm_rre(t, y, k):
    yprime = np.zeros_like(y)
    yprime[0] = -k[0] * y[0] * y[1] + k[1] * y[2]
    yprime[1] = -k[0] * y[0] * y[1] + (k[1] + k[2]) * y[2]
    yprime[2] = k[0] * y[0] * y[1] - (k[1] + k[2]) * y[2]
    yprime[3] = k[2] * y[2]

    return yprime

tspan = np.array([0.0, 50.0])
yzero = np.array([5.0e-7, 2.0e-7, 0.0, 0.0])
k = np.array([1.0e6, 1.0e-4, 0.1])
r = ode(mm_rre).set_integrator('vode', method='bdf', order=15)
r.set_f_params(k)
r.set_initial_value(yzero, tspan[0])
dt = 1.0
y = np.zeros((4,int(tspan[1]/dt)+1))
k = 1
while r.successful() and r.t < tspan[1]:
    r.integrate(r.t + dt)
    y[:,k] = r.y
    k+=1
    
tvals = np.linspace(tspan[0], tspan[1], int(tspan[1]/dt)+1)
plt.plot(tvals, y[0, :], 'go-', markersize=10, label='Substrate')
plt.plot(tvals, y[3, :], 'r*-', markersize=10, label='Product')
plt.xlabel("Time", size=14)
plt.ylabel("Concentration", size=14)
plt.legend(fontsize=12)
plt.show()

