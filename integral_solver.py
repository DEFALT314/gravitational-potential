import numpy as np
class GaussianQuadrature:
    def __init__(self, n):
        self.n = n
        self.nodes, self.weights = np.polynomial.legendre.leggauss(n)

    def integrate(self, func, a, b):
        integral = sum(((b - a) / 2) * self.weights[i] * func((b - a) / 2 * self.nodes[i] + (b + a) / 2) for i in range(self.n))
        return integral
    
    