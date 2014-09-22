import numpy as np
import matplotlib.pyplot as plt

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

count = 1
tvals = [0.0]
Xvals = [list(X)]

while t < tfinal:
    a[0] = c[0] * X[0] * X[1]
    a[1] = c[1] * X[2]
    a[2] = c[2] * X[2]
    asum = np.sum(a)
    j = np.argmax(np.random.rand() < np.cumsum(a / asum))
    tau = np.log(1.0 / np.random.rand()) / asum
    X += V[:, j]
    t += tau
    count += 1
    tvals.append(t)
    Xvals.append(list(X))

L = len(tvals)
tnew = np.zeros(2*L-1)
tnew[1:-1:2] = tvals[1:]
tnew[2::2] = tvals[1:]
tnew[0] = tvals[0]

Svals = np.array(Xvals)[:, 0]
ynew = np.zeros(2*L-1)
ynew[::2] = Svals
ynew[1:-1:2] = Svals[:-1]
plt.plot(tnew, ynew, 'go-', label = "Substrate")

Pvals = np.array(Xvals)[:, 3]
ynew = np.zeros(2*L-1)
ynew[::2] = Pvals
ynew[1:-1:2] = Pvals[:-1]
plt.plot(tnew, ynew, 'r*-', label = "Product")
plt.xlabel("Time", size=14)
plt.ylabel("Molecules", size=14)
plt.xlim((0.0, 55.0))
plt.ylim((0.0, 310.0))
plt.legend(fontsize=14)
plt.show()
