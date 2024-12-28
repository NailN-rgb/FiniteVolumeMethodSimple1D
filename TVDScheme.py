# implemenatation of TVD scheme for Convection-Diffusion Equation by QUICK Scheme
import numpy as np

# flux limiter function
def van_leer_limiter(r):
    return (r + abs(r)) / (1 + r)


def van_albada_limiter(r):
    return (r + pow(r,2)) / (1 + pow(r,2))


def quick_limiter(r):
    return max(0, min(2 * r, (3 + r) / 4, 2))


def min_mood_limiter(r):
    return min(r, 1) if r > 0 else 0


def superbee_limiter(r):
    return max(0, min(2 * r, 1), min(r, 2))


def sweby_limiter(r):
    beta = 1
    return max(0, min(beta * r, 1), min(r, beta))

def umsit_limiter(r):
    return max(
        0, 
        min(2 * r, (3 + r)/4, 2, (1 + 3 * r) / 4, 2)
    )


def init_phi(N):
    phi = []
    phinew = []
    for i in range(1, N  + 1):
        phi.append(i * 0.1)
        phinew.append(i * 0.1)

    return [phi, phinew]


def TVD_Scheme(N, x, equation, tol = 1e-4, tvd_limit = 10):
    A = np.zeros((N , N))
    F = np.zeros(N)

    limiter = umsit_limiter

    # phi arrays for iteration process
    [phi, phinew] = init_phi(N)

    for iteration in range(0, tvd_limit):
        # Fill matrix & rhs for nodes
        for i in range(1, N - 1):
            xE =  x[i + 1]
            xW =  x[i - 1]
            xP =  x[i]

            xe = (xP + xE) / 2
            xw = (xP + xW) / 2 

            Fe = equation.rho(xe) * equation.u(xe)
            Fw = equation.rho(xw) * equation.u(xw)

            De = equation.Gamma(xe) / (xE - xP)
            Dw = equation.Gamma(xw) / (xP - xW)

            alphae = (Fe > 0)
            alphaw = (Fw > 0)

            if i == 1:
                replus = (phi[i] - phi[i - 1]) / (phi[i + 1] - phi[i])
                rwplus = (phi[i - 1] - equation.ua) / (phi[i] - phi[i - 1])

                reminus = (phi[i + 2] - phi[i + 1]) / (phi[i + 1] - phi[i])
                rwminus = (phi[i + 1] - phi[i]) / (phi[i] - phi[i - 1])
            elif i == N - 2:
                replus = (phi[i] - phi[i - 1]) / (phi[i + 1] - phi[i])
                rwplus = (phi[i - 1] - phi[i - 2]) / (phi[i] - phi[i - 1])

                reminus = (equation.ub - phi[i + 1]) / (phi[i + 1] - phi[i])
                rwminus = (phi[i + 1] - phi[i]) / (phi[i] - phi[i - 1])
            else:
                replus = (phi[i] - phi[i - 1]) / (phi[i + 1] - phi[i])
                rwplus = (phi[i - 1] - phi[i - 2]) / (phi[i] - phi[i - 1])

                reminus = (phi[i + 2] - phi[i + 1]) / (phi[i + 1] - phi[i])
                rwminus = (phi[i + 1] - phi[i]) / (phi[i] - phi[i - 1])

            A[i, i - 1] = -(Dw + max(Fw, 0))
            A[i, i + 1] = -(De + max(-Fe, 0))
            A[i, i] = Dw + max(Fw, 0) + De + max(-Fe, 0) + (Fe - Fw)

            F[i] = 1/2 * Fe * ((1 - alphae) * limiter(reminus) - alphae * limiter(replus)) * (phi[i + 1] - phi[i]) + \
                            1/2 * Fw * (alphaw * limiter(rwplus) - (1 - alphaw) * limiter(rwminus)) * (phi[i] - phi[i - 1])

        # Boundary nodes
        # Node A
        D1 = equation.Gamma(x[1]) / (x[1] - x[0])
        Da = 2 * equation.Gamma(x[0]) / (x[1] - x[0])
        F1 = equation.rho(x[1]) * equation.u(x[1])
        Fa = equation.rho(x[0]) * equation.u(x[0])
        r1 = 2 * (phi[0] - equation.ua) / (phi[1] - phi[0])

        A[0, 0] = D1 + Da + F1
        A[0, 1] = -D1
        F[0]    = (Da + Fa) * equation.ua - 1/2 * Fe * limiter(r1) * (phi[1] - phi[0])

        # Node N
        Dn = equation.Gamma(x[-2]) / (x[-1] - x[-2])
        Db = 2 * equation.Gamma(x[-1]) / (x[-1] - x[-2])
        Fn = equation.rho(x[-2]) * equation.u(x[-2])
        Fb = equation.rho(x[-1]) * equation.u(x[-1])
        rn = 2 * (equation.ub - phi[-1]) / (phi[-1] - phi[-2])

        A[N - 1, N - 1] = Dn + Db
        A[N - 1, N - 2] = -Fn - Dn 
        F[N - 1]        = (Db - Fb) * equation.ub + 1/2 * Fn * limiter(rn) * (phi[-1] - phi[-2])

        phinew = np.linalg.solve(A,F)

        if(np.linalg.norm(phi - phinew) < tol):
            return phinew

        # print(A)
        # print(F)

        phi = phinew

        
    print("TVD-Iteration Limit")

    return phinew