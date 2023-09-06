import numpy as np

n = 10

a = -1.0
b = 2.0
c = -1.0

y = np.ones(n)

A = np.zeros((n, n))
A[0, 0] = b
A[0, 1] = c
for i in range(1, n - 1):
    A[i][i - 1] = a
    A[i][i] = b
    A[i][i + 1] = c

A[n - 1, n - 1] = b
A[n - 1, n - 2] = a

x = np.zeros_like(y)

for i in range(300):
    xtmp = x
    x[0] = (y[0] - c * xtmp[1]) / b
    for j in range(1,n-1):
        x[j] = (y[j] - a * xtmp[j-1] - c * xtmp[j+1]) / b
    x[n-1] = (y[n-1] - a * xtmp[n-2]) / b

    print(i,x)

print(np.dot(np.linalg.inv(A), y))