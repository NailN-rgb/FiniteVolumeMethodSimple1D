from steadydiffequatiom import SteadyStateEquation
from centraldiffcheme import central_diff_scheme
from steadydiffequatiom import get_modified_array


import matplotlib.pyplot as plt
import numpy as np

# determine equation
equation: SteadyStateEquation = SteadyStateEquation()


## create mesh
# Start & end point
a: float = 0.1
b: float = 2
# Nodes count
N: float = 20

x = np.linspace(a, b, N)

# get coorditanes of FV centers
x_mod = get_modified_array(x)

# Start solver
y = central_diff_scheme(N - 1, x, equation)

# Draw Solution plot
plt.plot(x_mod, y)
plt.savefig("solutionFVM.png")

