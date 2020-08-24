#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys


#====================================================
def make_2D_density_plot( xDir, yDir, xLimits, yLimits, \
                            xLabel, yLabel, nbins, logDensity ):
	fig, ax = plt.subplots()

	H, xedges, yedges = np.histogram2d(xDir, yDir, bins=nbins)

	#xcenters = (xedges[:-1] + xedges[1:]) / 2
	#ycenters = (yedges[:-1] + yedges[1:]) / 2

	cm = plt.cm.gnuplot
	im = None
	if logDensity:
		im = plt.imshow(H.T, cmap=cm, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], \
						norm=colors.LogNorm(vmin=(H[np.where(H>0.0)]).min(), vmax=H.max()), \
						origin='lower', interpolation='bicubic')
	else:
		im = plt.imshow(H.T, cmap=cm, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], \
						origin='lower', interpolation='bicubic')

	plt.colorbar(im,fraction=0.046, pad=0.04)

	plt.xlabel(xLabel)
	plt.ylabel(yLabel)

	#plt.axes().set_aspect('equal')
	plt.xlim(xLimits)
	plt.ylim(yLimits)
	
	plt.show(block = False)




#====================================================
def generate_plot( data, usingLogDensity, numberOfBins, \
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
		make_2D_density_plot(z, t, xLimits, yLimits, xLabel, yLabel, numberOfBins, usingLogDensity)
	elif plotMode == 1:
		make_2D_density_plot(x, y, xLimits, yLimits, xLabel, yLabel, numberOfBins, usingLogDensity)



#====================================================
if __name__ == "__main__":
	# Read in name of file from command line
	filename = sys.argv[1]
	
	# Load file
	data = np.loadtxt(filename).T
	
	# Generate plots
	generate_plot( data, False, 50, 1, [-5.0, 5.0], [-5.0, 5.0] )
	

# End of file
