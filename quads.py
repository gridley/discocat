#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

class EvenQuadrature:
    # Generates a level-symmetric S_n quadrature, given n/2
    # A nonlinear set of equations is created, with both
    # mu_1 and the weights of the quadrature being the unknowns

    def __init__(self, n_2):
        self.n_2 = n_2

    def residual(self, unknowns):
        # Unknowns are mu1, followed by the weights.
        mu1 = unknowns[0]
        weights = unknowns[1:]

        # Create residual
        residual = np.zeros(self.n_2+1)

        # Calculate mu values
        delta = 2.0 * (1.0 - 3.0 * mu1**2) / (2.0 * self.n_2 - 2.0)
        mus = np.zeros(self.n_2)
        mus[0] = mu1
        for i in range(1, self.n_2):
            mus[i] = np.sqrt(mus[i-1]**2 + delta)

        print(mus)

        for i in range(self.n_2+1):
            residual[i] = 2.0 * np.dot(weights, mus**(2*i)) - 1.0/(2.0 * i + 1)

        return residual

e4 = EvenQuadrature(2)
print(e4.residual([0.3500212, 0.3333333, 0.1666667]))
