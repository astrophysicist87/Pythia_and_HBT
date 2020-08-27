#!/usr/bin/env python
#====================

import numpy as np
import sys

#====================================================
if __name__ == "__main__":
    # Read in name of file from command line
    filename = sys.argv[1]
    KT = float(sys.argv[2])
	
    # Load file
    data = np.loadtxt(filename).T
    [Kphi,t,x,y,z]=data
    xo = x*np.cos(Kphi) + y*np.sin(Kphi)
    xs = -x*np.sin(Kphi) + y*np.cos(Kphi)
    xl = z
    
    mpi = 0.13957
    betaT = KT / np.sqrt(mpi**2 + KT**2)
    
    print '<x_o>, <x_s>, <x_l>, <t>:'
    print np.mean(xo), np.mean(xs), np.mean(xl), np.mean(t)
    print '<(x_o-<x_o>)^2>, <(x_s-<x_s>)^2>, <(x_l-<x_l>)^2>, <(x_o-<x_o>)(t-<t>)>, <(t-<t>)^2>:'
    print np.var(xo), np.var(xs), np.var(xl), \
            np.mean(xo*t)-np.mean(xo)*np.mean(t), np.var(t)
    print 'R^2_o, R^2_s, R^2_l:'
    print np.var(xo) - 2.0*betaT*(np.mean(xo*t)-np.mean(xo)*np.mean(t))+betaT**2*np.var(t),\
          np.var(xs), np.var(xl)


# End of file
