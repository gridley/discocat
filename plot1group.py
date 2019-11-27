#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
dat = np.loadtxt('flux')
l = dat.shape[0]
mesh_dimx = int(l**0.5)
dat2 = dat.reshape((mesh_dimx, mesh_dimx))
plt.pcolormesh(dat2)
plt.colorbar()
plt.savefig('1group_%i.png' % mesh_dimx)
