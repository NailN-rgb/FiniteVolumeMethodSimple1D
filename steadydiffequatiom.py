import numpy as np
# For equation 
# dib(Gamma (grad phi)) + S_phi = 0

class SteadyStateEquation:
    # 'Flow' 
    def Gamma(self, x):
        return x
    
    # 'Clean' Source
    def Su(self, x):
        return 0
    
    # 'Node' Source
    def Sp(self, x):
        return 0
    
    # density
    def rho(self, x):
        return 0
    
    # flow field
    def u(self, x):
        return 0
    
    ua = -2.30
    ub = 0.69


def get_modified_array(x):
    result = []

    for i in range(1, x.size):
        result.append((x[i] + x[i-1]) / 2)

    return result