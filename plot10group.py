#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
dat = np.loadtxt('flux')
l = dat.shape[0] / 10
mesh_dimx = int(l**0.5)
for g in range(10):
    dat2 = dat[g::10].reshape((mesh_dimx, mesh_dimx))
    plt.pcolormesh(dat2)
    plt.title('Group %i' % g)
    plt.colorbar()
    plt.show()
