import numpy as np
# For equation 
# dib(Gamma (grad phi)) + S_phi = 0

class SteadyStateEquation:
    def Gamma(self, x):
        return x
    
    def Su(self, x):
        return 0
    
    def Sp(self, x):
        return 0
    
    ua = -2.30
    ub = 0.69


def get_modified_array(x):
    result = []

    for i in range(1, x.size):
        result.append((x[i] + x[i-1]) / 2)

    return result