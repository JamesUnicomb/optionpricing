import numpy as np

n = 31
a = -1.0 * np.ones(n)
a[0] = 0.0  # fix indexing
b = 2.0 * np.ones(n)
c = -1.0 * np.ones(n)
c[n - 1] = 0.0  # fix indexing

y = np.ones(n)

A = np.zeros((n, n))
A[0, 0] = b[0]
A[0, 1] = c[0]
for i in range(1, n - 1):
    A[i][i - 1] = a[i]
    A[i][i] = b[i]
    A[i][i + 1] = c[i]

A[n - 1, n - 1] = b[n - 1]
A[n - 1, n - 2] = a[n - 2]

print("A = ", A)
print("y = ", y)
print("x = ", np.dot(np.linalg.inv(A), y))

# start the algorithm
cprime = np.zeros_like(c)

cprime[0] = c[0] / b[0]
for i in range(1, len(c)):
    cprime[i] = c[i] / (b[i] - cprime[i - 1] * a[i])

yprime = np.zeros_like(y)

yprime[0] = y[0] / b[0]
for i in range(1, len(y)):
    yprime[i] = (y[i] - yprime[i - 1] * a[i]) / (b[i] - cprime[i - 1] * a[i])

xprime = np.zeros(n)

xprime[n - 1] = yprime[n - 1]

for i in range(n - 2, -1, -1):
    xprime[i] = yprime[i] - cprime[i] * xprime[i + 1]

print("xprime = ", xprime)
