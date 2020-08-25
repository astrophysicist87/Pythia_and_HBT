#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

#====================================================
def pause():
    programPause = raw_input()

#====================================================
def make_2D_density_plot( xDir, yDir, xLimits, yLimits, \
                            xLabel, yLabel, bws ):
    [xbinwidth, ybinwidth] = bws
    nxbins = int((xLimits[1] - xLimits[0]) / xbinwidth)
    nybins = int((yLimits[1] - yLimits[0]) / ybinwidth)
    fig, ax = plt.subplots()

    try:
        H, xedges, yedges = np.histogram2d(xDir, yDir, bins=[nxbins, nybins],
                                        range=[xLimits, yLimits], density=True)
    except:
        H, xedges, yedges = np.histogram2d(xDir, yDir, bins=[nxbins, nybins],
                                        range=[xLimits, yLimits], normed=True)

    cm = plt.cm.gnuplot
    im = plt.imshow(H.T, cmap=cm, origin='lower', interpolation='bilinear',
                    extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()])

    plt.colorbar(im,fraction=0.046, pad=0.04)

    plt.xlabel(xLabel)
    plt.ylabel(yLabel)

    plt.xlim(xLimits)
    plt.ylim(yLimits)
	
    plt.show(block = False)




#====================================================
def generate_plot( data, polarMode, bws, \
                    plotMode, xLimits, yLimits ):
    if plotMode == 0:
        xLabel = r'$z$ (fm)'
        yLabel = r'$t$ (fm/$c$)'
        if polarMode == 1:
            xLabel = r'$\eta$'
            yLabel = r'$\tau$ (fm/$c$)'
    elif plotMode == 1:
        xLabel = r'$x$ (fm)'
        yLabel = r'$y$ (fm)'
        if polarMode == 1:
            xLabel = r'$r$ (fm)'
            yLabel = r'$\phi$'
    elif plotMode == 2:
        xLabel = r'$r$ (fm)'
        yLabel = r'$\tau$ (fm)'
    else:
        print 1/0

    [Kphi,t,x,y,z]=data
	
    # Plot them
    if plotMode == 0:
        if polarMode == 0:
            make_2D_density_plot(z, t, xLimits, yLimits, xLabel, yLabel, bws)
        else:
            safeIndices = np.where(np.greater(t**2, z**2))
            tsI = t[safeIndices]
            zsI = z[safeIndices]
            eta = 0.5*np.log( (tsI+zsI)/(tsI-zsI) )
            tau = np.sqrt(tsI**2-zsI**2)
            make_2D_density_plot(eta, tau, xLimits, yLimits, xLabel, yLabel, bws)
    elif plotMode == 1:
        make_2D_density_plot(x, y, xLimits, yLimits, xLabel, yLabel, bws)
    elif plotMode == 2:
        safeIndices = np.where(np.greater(t**2, z**2))
        tsI = t[safeIndices]
        xsI = x[safeIndices]
        ysI = y[safeIndices]
        zsI = z[safeIndices]
        r = np.sqrt(xsI**2+ysI**2)
        tau = np.sqrt(tsI**2-zsI**2)
        make_2D_density_plot(r, tau, xLimits, yLimits, xLabel, yLabel, bws)



#====================================================
if __name__ == "__main__":
    # Read in name of file from command line
    filename = sys.argv[1]
	
    # Load file
    data = np.loadtxt(filename).T
	
    # Generate plots
    generate_plot( data, 0, [0.1, 0.1], 0, [-10.0, 10.0], [0.0, 15.0] )
    generate_plot( data, 1, [0.1, 0.1], 0, [-5.0, 5.0], [0.0, 5.0] )
    generate_plot( data, 0, [0.1, 0.1], 1, [-5.0, 5.0], [-5.0, 5.0] )
    generate_plot( data, None, [0.1, 0.1], 2, [0.0, 5.0], [0.0, 5.0] )


    pause()

# End of file