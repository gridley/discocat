#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
dat = np.loadtxt(sys.argv[1])
l = dat.shape[0] / 10
mesh_dimx = int(l**0.5)
for g in range(10):
    dat2 = dat[g::10].reshape((mesh_dimx, mesh_dimx))
    plt.pcolormesh(dat2)
    plt.title('Group %i' % g)
    if g == 0:
        plt.colorbar()
    plt.savefig('group%i.png' % g)
