import numpy as np

# Stoichiometric matrix
V = np.array([[-1.0, 1.0, 0.0],[-1.0, 1.0, 1.0],[1.0, -1.0, -1.0],[0.0, 0.0, 1.0]])

# Parameters and Initial Conditions
nA = 6.023e23             # Avagadro's number
vol = 1e-15               # volume of system
Y = np.zeros((4,))
c = np.zeros((3,))
d = np.zeros_like(c)
a = np.zeros_like(c)
Y[0] = round(5e-7 * nA * vol) # molecules of substrate
Y[1] = round(2e-7 * nA * vol) # molecules of enzyme 
c[0] = 1.0e6 / (nA * vol)
c[1] = 1.0e-4
c[2] = 0.1

tfinal = 50.0
L = 250
tau = tfinal / L    # stepsize

for k in range(L):
    a[0] = c[0] * Y[0] * Y[1]
    a[1] = c[1] * Y[2]
    a[2] = c[2] * Y[2]
    for i in range(3):
        d[i] = tau * a[i] + np.sqrt(abs(tau * a[i])) * np.random.randn()
    Y = Y + d[0] * V[:, 0] + d[1] * V[:, 1] + d[2] * V[:, 2]
