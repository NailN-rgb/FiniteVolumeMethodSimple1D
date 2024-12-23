import numpy as np

# Central differencing scheme 
# Applied for 1D Diffusion problems

def central_diff_scheme(N, x, equation):
    A = np.zeros((N , N))
    F = np.zeros(N)

    # Fill matrix for nodes
    for i in range(1, N - 1):
        xe = x[i+1]
        xw = x[i]
        A[i,i]    =  equation.Gamma(xe) / (xe - xw) + equation.Gamma(xw) / (xe - xw)  - equation.Sp((xe + xw) / 2)
        A[i, i-1] = -equation.Gamma(xw) / (xe - xw) 
        A[i, i+1] = -equation.Gamma(xe) / (xe - xw)

        F[i] = equation.Su((xe + xw) / 2)

    # Assemble first BC for left end
    A[0, 0] = (2 * equation.Gamma(x[0]) + equation.Gamma(x[1])) / (x[1] - x[0])
    A[0, 1] =     -equation.Gamma(x[1]) / (x[1] - x[0])
    F[0]    =  2 * equation.Gamma(x[0]) / (x[1] - x[0]) * equation.ua 

    # Assemble first BC for right end
    A[N - 1, N - 1] = (2 * equation.Gamma(x[-1]) + equation.Gamma(x[-2])) / (x[-1] - x[-2])
    A[N - 1, N - 2] =     -equation.Gamma(x[-2]) / (x[-1] - x[-2])
    F[N - 1]        =  2 * equation.Gamma(x[-1]) / (x[-1] - x[-2]) * equation.ub

    y = np.linalg.solve(A,F)
    return y