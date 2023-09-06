import numpy as np
from timeit import timeit
import optionpricing

n = 32767
a = -1.0
b = 2.0
c = -1.0

solver = optionpricing.TriDiagSolver(n, a, b, c)

y = [1.0 for i in range(n)]

# solver.sety(y)
# solver.solve()
# #print(np.round(solver.getx(), 3).tolist())

# solver.sety(y)
# solver.solve()
# #print(np.round(solver.getx(), 3).tolist())


def f1():
    solver.sety(y)
    solver.solve_parallel()

def f2():
    solver.sety(y)
    solver.solve_serial()


print("parallel: ", timeit(f1, number=1000) / 1000)
print("serial: ", timeit(f2, number=1000) / 1000)

# y = np.ones(n)

# A = np.zeros((n, n))
# A[0, 0] = b
# A[0, 1] = c
# for i in range(1, n - 1):
#     A[i][i - 1] = a
#     A[i][i] = b
#     A[i][i + 1] = c

# A[n - 1, n - 1] = b
# A[n - 1, n - 2] = a

# #print(np.round(np.dot(np.linalg.inv(A), y), 3).tolist())