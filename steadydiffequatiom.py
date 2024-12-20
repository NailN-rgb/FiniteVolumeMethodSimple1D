
# For equation 
# dib(Gamma (grad phi)) + S_phi = 0

class SteadyStateEquation:
    Gamma = 1
    Su = 0
    Sp = 0
    ua = 0
    ub = 1


def get_modified_array(x):
    result = []

    for i in range(1, x.size):
        result.append((x[i] + x[i-1]) / 2)

    return result