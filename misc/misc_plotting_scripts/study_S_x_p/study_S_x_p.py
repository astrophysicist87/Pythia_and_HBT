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
    xLabel, yLabel, nbins ):
    fig, ax = plt.subplots()

    H, xedges, yedges = np.histogram2d(xDir, yDir, bins=nbins, range=[xLimits, yLimits])

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
def generate_plot( data, numberOfBins, \
                    plotMode, xLimits, yLimits ):
    if plotMode == 0:
        xLabel = r'$z$ (fm)'
        yLabel = r'$t$ (fm/$c$)'
    elif plotMode == 1:
        xLabel = r'$x$ (fm)'
        yLabel = r'$y$ (fm)'
    else:
        print 1/0

    [t,x,y,z]=data
	
    # Plot them
    if plotMode == 0:
        make_2D_density_plot(z, t, xLimits, yLimits, xLabel, yLabel, numberOfBins)
    elif plotMode == 1:
        make_2D_density_plot(x, y, xLimits, yLimits, xLabel, yLabel, numberOfBins)



#====================================================
if __name__ == "__main__":
    # Read in name of file from command line
    filename = sys.argv[1]
	
    # Load file
    data = np.loadtxt(filename).T
	
    # Generate plots
    generate_plot( data, 50, 0, [-10.0, 10.0], [0.0, 15.0] )
    generate_plot( data, 50, 1, [-5.0, 5.0], [-5.0, 5.0] )


    pause()

# End of file
