#!/usr/bin/env python3
# How do we know that the solution to these nonlinear equations is unique?
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as lin
import scipy.optimize
from mpl_toolkits.mplot3d import Axes3D
import itertools

class LevelSymmetricQuadrature:

    # This class defines symmetry relations required for pointwise weights
    # on a quadrature set, and provides a plotting utility.

    # Level weight to point weight matrices (these may be singular!)
    symmetry_matrices = {
    2:      [[1.0]],

    4:      [[2.0, 0.0],
             [1.0, 0.0]],

    6:      [[2.0, 1.0],
             [0.0, 2.0],
             [1.0, 0.0]],

    8:      [[2.0, 2.0, 0.0],
             [0.0, 2.0, 1.0],
             [0.0, 2.0, 0.0],
             [1.0, 0.0, 0.0]],

    10:      [[2.0, 2.0, 1.0, 0.0],
              [0.0, 2.0, 0.0, 2.0],
              [0.0, 0.0, 2.0, 1.0],
              [0.0, 2.0, 0.0, 0.0],
              [1.0, 0.0, 0.0, 0.0]],

    12:      [[2.0, 2.0, 2.0, 0.0, 0.0],
              [0.0, 2.0, 0.0, 2.0, 1.0],
              [0.0, 0.0, 2.0, 0.0, 2.0],
              [0.0, 0.0, 2.0, 1.0, 0.0],
              [0.0, 2.0, 0.0, 0.0, 0.0],
              [1.0, 0.0, 0.0, 0.0, 0.0]],

    14:      [[2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0],
              [0.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0],
              [0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 1.0],
              [0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 0.0],
              [0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0],
              [0.0, 0.0,-1.0, 1.0, 1.0,-2.0, 1.0],
              [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],

    16:       [[2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0],
              [0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 2.0],
              [0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 1.0],
              [0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0],
              [0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0],
              [0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
              #[0.0, 0.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0]]
    }

    def plot(self):
        # Create a 3D plot of the quadrature
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Loop over all points on the sphere (in positive quadrant)
        n = self.n_2 * 2
        npts = int(n*(n+2)/8) # this is guaranteed to be an integer
        xvals = np.zeros(npts)
        yvals = np.zeros(npts)
        zvals = np.zeros(npts)
        weights = np.zeros(npts) 

        pt = 0
        ptweight_dict = {}
        symm_group = 0
        for i in range(self.n_2):
            for j in range(self.n_2-i):
                k = self.n_2-i-j-1

                # Update dictionary of point weights
                # Each different point weight corresponds to a group
                # of permuations of the indices. This comes from the
                # symmetry condition on point weights.
                if (i, j, k) not in ptweight_dict:
                    newkeys = itertools.permutations((i,j,k))
                    for newkey in newkeys:
                        ptweight_dict[newkey] = symm_group
                    symm_group += 1

                xvals[pt] = self.mus[i]
                yvals[pt] = self.mus[j]
                zvals[pt] = self.mus[k]
                weights[pt] = self.point_weights[ptweight_dict[(i,j,k)]]
                pt += 1

        # permute scatter plot over all octants
        # The below for loops can be commented and the next line uncommented
        # to get a single octant.
        # xs = 1; ys = 1; zs = 1
        sgns = [-1.0, 1.0]
        for xs in sgns:
            for ys in sgns:
                for zs in sgns:
                    sc = ax.scatter(xs*xvals, ys*yvals, zs=zs*zvals, s=20, c=weights)
        plt.colorbar(sc)

        # Plot some circles that intersect the points in order to exhibit the
        # unique nature of the level symmetric quadratures:
        n_circ = 50 # points on the circle
        phis = np.linspace(0.0, 2.0 * np.pi, n_circ)
        for i in range(self.n_2):
            xvals = np.zeros(n_circ)
            yvals = np.zeros(n_circ)
            zvals = np.zeros(n_circ)
            j = 0
            k = self.n_2-i-j-1
            circ_radius = np.sqrt(self.mus[j]**2 + self.mus[k]**2)
            for l in range(n_circ):
                xvals[l] = np.abs(self.mus[i])
                yvals[l] = circ_radius * np.cos(phis[l])
                zvals[l] = circ_radius * np.sin(phis[l])
            for sign in [-1, 1]:
                ax.plot(sign * xvals, yvals, zvals, c='b', alpha=0.25)
            for sign in [-1, 1]:
                ax.plot(yvals, sign*xvals, zvals, c='b', alpha=0.25)
            for sign in [-1, 1]:
                ax.plot(yvals, zvals, sign*xvals, c='b', alpha=0.25)

        plt.show()

    def calcPointWeights(self):
        # assumes that self.mus and self.weights have already been calculated
        n = self.n_2 * 2
        
        # Division by two here matches the normalization of the paper
        mat_sym = np.array(LevelSymmetricQuadrature.symmetry_matrices[n]) / 2.0
        self.point_weights = np.matmul(lin.pinv(mat_sym),self.weights)

class EvenQuadrature(LevelSymmetricQuadrature):
    # Generates a level-symmetric S_n quadrature, given n
    # A nonlinear set of equations is created, with both
    # mu_1 and the weights of the quadrature being the unknowns
    # This quadrature set exactly integrates even power monomials 
    # of the polar angle.
    #
    # This then feeds into a linear system for the weights,
    # which arises from symmetry considerations on levels.

    def __init__(self, n):
        self.n_2 = int(n/2)

        # avoid writing "self" a bunch
        n_2 = self.n_2

        # solve for quadrature
        guess = np.zeros(n_2 + 1)
        guess[0] = 0.10 # first cosine
        guess[1:] = np.ones(n_2) / n_2 # equal weights
        soln = scipy.optimize.root(self.residual, guess)

        if not soln.success:
            raise Exception('Failed to find quadrature in nonlinear solve')

        self.mus = np.zeros(self.n_2)
        self.mus[0] = soln.x[0]
        delta = 2.0 * (1.0 - 3.0 * self.mus[0]**2) / (2.0 * self.n_2 - 2.0)
        for i in range(1, self.n_2):
            self.mus[i] = np.sqrt(self.mus[i-1]**2 + delta)
        self.weights = soln.x[1:]

        self.calcPointWeights()

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

        for i in range(self.n_2+1):
            residual[i] = 2.0 * np.dot(weights, mus**(2*i)) - 1.0/(2.0 * i + 1)

        return residual

e4 = EvenQuadrature(12)
print(e4.mus)
print(e4.weights)
print(e4.point_weights)
e4.plot()
