import numpy as np

def central_diff_scheme(N, x, equation):
    A = np.zeros((N - 1, N - 1))
    F = np.zeros(N - 1)

    # Fill matrix for nodes
    for i in range(1, N-2):
        A[i,i]    =  (equation.Gamma) / (x[i + 1] - x[i]) + (equation.Gamma) / (x[i] - x[i-i]) - equation.Sp
        A[i, i-1] = -(equation.Gamma) / (x[i + 1] - x[i]) 
        A[i, i+1] = -(equation.Gamma) / (x[i] - x[i - 1])

        F[i] = equation.Su

    # Assemble first BC
    A[0, 0] = 2 * (equation.Gamma) / (x[i + 1] - x[i])
    F[0]    = equation.ua

    A[N - 2, N - 2] = 2 * (equation.Gamma) / (x[i + 1] - x[i])
    F[N - 2]        = equation.ub

    y = np.linalg.solve(A,F)
    return y