from steadydiffequatiom import SteadyStateEquation
from centraldiffcheme import central_diff_scheme
from steadydiffequatiom import get_modified_array


import matplotlib.pyplot as plt
import numpy as np

# determine equation
equation: SteadyStateEquation = SteadyStateEquation()


## create mesh
# Start & end point
a: float = 0
b: float = 1
# Nodes count
N: float = 10

x = np.linspace(a, b, N)
x_mod = get_modified_array(x)

y = central_diff_scheme(N, x, equation)
# Assemble matrix
plt.plot(x, y)
plt.savefig("solutionFVM.png")

