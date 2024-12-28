from steadydiffequatiom import SteadyStateEquation
from centraldiffcheme import central_diff_scheme
from steadydiffequatiom import get_modified_array
from TVDScheme import TVD_Scheme


import matplotlib.pyplot as plt
import numpy as np

# determine equation
equation: SteadyStateEquation = SteadyStateEquation()


## create mesh
# Start & end point
a: float = 0
b: float = 1
# Nodes count
N: float = 100

x = np.linspace(a, b, N)

# get coorditanes of FV centers
x_mod = get_modified_array(x)

# Start solver
y = TVD_Scheme(N - 1, x, equation)
y_cor = equation.analytical_solution(x_mod)
# Draw Solution plot
plt.plot(x_mod, y)
plt.plot(x_mod, y_cor, 'r')
plt.savefig("solutionFVM.png")

