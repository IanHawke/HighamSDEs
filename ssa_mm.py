import numpy as np

# Stoichiometric matrix
V = np.array([[-1.0, 1.0, 0.0],[-1.0, 1.0, 1.0],[1.0, -1.0, -1.0],[0.0, 0.0, 1.0]])

# Parameters and Initial Conditions
nA = 6.023e23             # Avagadro's number
vol = 1e-15               # volume of system
X = np.zeros((4,))
c = np.zeros((3,))
d = np.zeros_like(c)
a = np.zeros_like(c)
X[0] = round(5e-7 * nA * vol) # molecules of substrate
X[1] = round(2e-7 * nA * vol) # molecules of enzyme 
c[0] = 1.0e6 / (nA * vol)
c[1] = 1.0e-4
c[2] = 0.1

t = 0.0
tfinal = 50.0

while t < tfinal:
    a[0] = c[0] * X[0] * X[1]
    a[1] = c[1] * X[2]
    a[2] = c[2] * X[2]
    asum = np.sum(a)
    j = np.argmax(np.random.rand() < np.cumsum(a / asum))
    tau = np.log(1.0 / np.random.rand()) / asum
    X += V[:, j]
    t += tau
