import numpy as np

import optionpricing


def test_solve():
    n = 8
    a = 1.0
    b = 2.0
    c = 1.0

    A = np.zeros((n, n))
    A[0, 0] = b
    A[0, 1] = c
    for i in range(1, n - 1):
        A[i][i - 1] = a
        A[i][i] = b
        A[i][i + 1] = c
    A[n - 1, n - 1] = b
    A[n - 1, n - 2] = a

    print(A)

    y = np.arange(n)

    print(np.dot(np.linalg.inv(A), y))

    solver = optionpricing.TriDiagSolver(n, a, b, c)
    solver.sety(y)

    solver.solve()

    print(np.linalg.inv(A))
    print(solver.getx())
