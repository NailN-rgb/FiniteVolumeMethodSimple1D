import numpy as np
# For equation 
# d/dx (pu\phi) = d/dx(Г dphi/dx)
class SteadyStateEquation:
    # 'Flow' 
    def analytical_solution(self, x):
        res = []
        for point in x:
            res.append(-(np.exp((self.rho(point) * self.u(point) * point) / self.Gamma(point)) - 1) / \
                        (np.exp((self.rho(point) * self.u(point))/self.Gamma(point)) - 1) + 1)
        
        return res

    def Gamma(self, x):
        return 0.1
    
    # 'Clean' Source
    def Su(self, x):
        return 0
    
    # 'Node' Source
    def Sp(self, x):
        return 0
    
    # density
    def rho(self, x):
        return 1
    
    # flow field
    def u(self, x):
        return -2.5
    
    ua = 1
    ub = 0


def get_modified_array(x):
    result = []

    for i in range(1, x.size):
        result.append((x[i] + x[i-1]) / 2)

    return result