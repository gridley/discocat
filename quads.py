#!/usr/bin/env python3
# This generates level-symmetric quadratures for an S_n solver
# that match those obtained by Lathrop and Carlson.
#
# Note:
# How do we know that the solution to these nonlinear equations is unique?
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as lin
import scipy.optimize
from mpl_toolkits.mplot3d import Axes3D
import itertools
from textwrap import wrap

np.seterr(all='raise')

# helper to wrap in C++ initializer syntax
def wrapCurly(ll):
    lls = [str(s) for s in ll]
    return '\n      '.join(wrap('{' + ', '.join(lls) + '}'))

class LevelSymmetricQuadrature:
    '''
    This class defines symmetry relations required for pointwise weights
    on a quadrature set, and provides a plotting utility.
    '''


    def generate_symmetry_dict(self):
        '''
        Generates a dictionary that maps i,j,k indices of the
        quadrature to its enumerated group that tells which
        rays must have the same weight in order to maintain
        symmetry about all 90 deg rotations
        '''

        if self.ptweight_dict:
            return

        symm_group = 0
        for i in range(self.n_2):
            for j in range(self.n_2-i):
                k = self.n_2-i-j-1

                # Update dictionary of point weights
                # Each different point weight corresponds to a group
                # of permuations of the indices. This comes from the
                # symmetry condition on point weights.
                if (i, j, k) not in self.ptweight_dict:
                    newkeys = itertools.permutations((i,j,k))
                    for newkey in newkeys:
                        self.ptweight_dict[newkey] = symm_group
                    self.max_symm_group = symm_group
                    symm_group += 1

    def gen_sym(self):
        '''
        Generates a singular matrix, which, if the RHS satisfied certain conditions, 
        yields a system of equations where level weights are the RHS, and
        the solution of this system gives a set of weights for each ray in
        a symmetric group.

        e.g.
                         1
           1            2 2
          2 2          3 4 3
         2 3 2        2 4 4 2
        1 2 2 1      1 2 3 2 1

        are each equivalent to a certain matrix
        '''

        self.generate_symmetry_dict()

        ncols = max(self.ptweight_dict.values()) + 1
        nrows = self.n_2

        # Allocate matrix:
        mat = np.zeros((nrows, ncols))

        self.generate_symmetry_dict()
        for i in range(self.n_2):
            for j in range(self.n_2-i):
                k = self.n_2-i-j-1
                group = self.ptweight_dict[(i,j,k)]
                mat[i, group] += 1.0
        return mat

    def plot(self, circles=False):

        # Make sure symmetry dictionary is created:
        self.generate_symmetry_dict()

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
        for i in range(self.n_2):
            for j in range(self.n_2-i):
                k = self.n_2-i-j-1
                xvals[pt] = self.mus[i]
                yvals[pt] = self.mus[j]
                zvals[pt] = self.mus[k]
                weights[pt] = self.point_weights[self.ptweight_dict[(i,j,k)]]
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
        if circles:
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
    def calcPointWeights(self, sum4pi=False):
        # assumes that self.mus and self.weights have already been calculated
        n = self.n_2 * 2
        
        # Division by two here matches the normalization of the paper
        mat_sym = self.gen_sym() / 2.0
        self.point_weights = np.matmul(lin.pinv(mat_sym),self.weights)

        if sum4pi:
            # Get count of 
            tot_weight = np.sum(self.point_weights)
            self.point_weights *= 4.0 * np.pi / tot_weight / (n * (n+2))

    def toCpp(self):
        # Prints out a template specialization of my LevelSymmetricQuadrature class
        # that gives a nice, clean syntax.

        retstr = """
  else if (na == %i and ma == %s)
  {
    mu = %s;
    weights = %s;
    point_weights = %s;
  }
""" % ( self.n_2 * 2,
        self.type,
        wrapCurly(self.mus),
        wrapCurly(self.weights),
        wrapCurly(self.point_weights))
        return retstr

class EvenQuadrature(LevelSymmetricQuadrature):
    # Generates a level-symmetric S_n quadrature, given n
    # A nonlinear set of equations is created, with both
    # mu_1 and the weights of the quadrature being the unknowns
    # This quadrature set exactly integrates even power monomials 
    # of the polar angle.
    #
    # This then feeds into a linear system for the weights,
    # which arises from symmetry considerations on levels.

    def __init__(self, n, sum4pi=False):
        self.type = 'EVEN'
        self.n_2 = int(n/2)

        # avoid writing "self" a bunch
        n_2 = self.n_2

        # empty point weight group dictionary
        self.ptweight_dict = {}

        # solve for quadrature
        guess = np.zeros(n_2 + 1)
        guess[1:] = np.ones(n_2) / n_2 # equal weights

        
        while True:
            try:
                guess[0] = 0.5 * np.random.rand(1)[0]
                soln = scipy.optimize.root(self.residual, guess, method='anderson')
                if not soln.success:
                    continue
                break
            except:
                pass

        if not soln.success:
            raise Exception('Failed to find quadrature in nonlinear solve')

        self.mus = np.zeros(self.n_2)
        self.mus[0] = np.abs(soln.x[0])
        delta = 2.0 * (1.0 - 3.0 * self.mus[0]**2) / (2.0 * self.n_2 - 2.0)
        for i in range(1, self.n_2):
            self.mus[i] = np.sqrt(self.mus[i-1]**2 + delta)
        self.weights = soln.x[1:]

        self.calcPointWeights(sum4pi=sum4pi)

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

class OddQuadrature(LevelSymmetricQuadrature):
    # Generates a level-symmetric S_n quadrature, given n
    # a quadrature which matches odd moments.

    def __init__(self, n, sum4pi=False):
        self.type = 'ODD'
        self.n_2 = int(n/2)

        # avoid writing "self" a bunch
        n_2 = self.n_2

        # empty point weight group dictionary
        self.ptweight_dict = {}

        # solve for quadrature
        guess = np.zeros(n_2 + 1)
        guess[1:] = np.ones(n_2) / n_2 # equal weights

        while True:
            try:
                guess[0] = 0.5 * np.random.rand(1)[0]
                soln = scipy.optimize.root(self.residual, guess, method='anderson')
                break
            except:
                pass

        if not soln.success:
            raise Exception('Failed to find quadrature in nonlinear solve')

        self.mus = np.zeros(self.n_2)
        self.mus[0] = np.abs(soln.x[0])
        delta = 2.0 * (1.0 - 3.0 * self.mus[0]**2) / (2.0 * self.n_2 - 2.0)
        for i in range(1, self.n_2):
            self.mus[i] = np.sqrt(self.mus[i-1]**2 + delta)
        self.weights = soln.x[1:]

        self.calcPointWeights(sum4pi=sum4pi)

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
            residual[i] = 2.0 * np.dot(weights, mus**i) - 1.0/(i + 1)

        return residual

# for i in range(4, 30, 2):
#     even = EvenQuadrature(i)#, sum4pi=True)
#     odd = OddQuadrature(i)#, sum4pi=True)
#     print(even.toCpp())
#     print(odd.toCpp())

# print(q.gen_sym())
# print(q.gen_sym())
q = EvenQuadrature(4, sum4pi=True)
# print(q.toCpp())
# print(q.toCpp())
q.plot()
