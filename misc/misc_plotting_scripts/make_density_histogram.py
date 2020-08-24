#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

# Useful constants, settings, etc.
mmTOfm          = 1.e12

usingLogDensity = False
numberOfBins    = 100

# Plotting settings ( x and y refer to figure axes,
#					  not necessarily physical coordinates )
xLimits = [-10.0,10.0]
yLimits = [-10.0,10.0]


#====================================================
def pause():
    programPause = raw_input()


#====================================================
def make_2D_density_plot_v1( xDir, yDir, xLabel='', yLabel='' ):
	fig, ax = plt.subplots()

	H, xedges, yedges = np.histogram2d(xDir, yDir, bins=numberOfBins)

	xcenters = (xedges[:-1] + xedges[1:]) / 2
	ycenters = (yedges[:-1] + yedges[1:]) / 2

	cm = plt.cm.gnuplot
	im = None
	if usingLogDensity:
		im = plt.imshow(H.T, cmap=cm, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], \
						norm=colors.LogNorm(vmin=(H[np.where(H>0.0)]).min(), vmax=H.max()), \
						origin='lower', interpolation='bicubic')
	else:
		im = plt.imshow(H.T, cmap=cm, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], \
						origin='lower', interpolation='bicubic')

	plt.colorbar(im,fraction=0.046, pad=0.04)
	ax.patch.set_facecolor('black')

	plt.xlabel(xLabel)
	plt.ylabel(yLabel)

	#plt.axes().set_aspect('equal')
	plt.xlim(xLimits)
	plt.ylim(yLimits)
	
	#plt.show(block = False)
	plt.savefig("./posICs.pdf", bbox_inches='tight')

#====================================================
def make_2D_density_plot_v2( xDir, yDir, zDir, xLabel='', yLabel='' ):
	fig, ax = plt.subplots()

	H, xedges, yedges = np.histogram2d(xDir, yDir, weights=zDir, bins=numberOfBins)

	xcenters = (xedges[:-1] + xedges[1:]) / 2
	ycenters = (yedges[:-1] + yedges[1:]) / 2

	cm = plt.cm.gnuplot
	im = None
	if usingLogDensity:
		im = plt.imshow(H.T, cmap=cm, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], \
						norm=colors.LogNorm(vmin=(H[np.where(H>0.0)]).min(), vmax=H.max()), \
						origin='lower', interpolation='bicubic')
	else:
		im = plt.imshow(H.T, cmap=cm, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], \
						origin='lower', interpolation='bicubic')

	plt.colorbar(im,fraction=0.046, pad=0.04)
	ax.patch.set_facecolor('black')

	plt.xlabel(xLabel)
	plt.ylabel(yLabel)

	#plt.axes().set_aspect('equal')
	plt.xlim(xLimits)
	plt.ylim(yLimits)
	
	#plt.show(block = False)
	plt.savefig("./ewICs.pdf", bbox_inches='tight')


#====================================================
if __name__ == "__main__":
	# Read in name of file and columns to plot from command line
	filename = sys.argv[1]
	cols = sys.argv[2:]
	
	# Load file
	data = np.loadtxt(filename, usecols=tuple(map(lambda z: int(z), cols)))
	#data = data[ np.where( (np.abs(data[:,0]) > 0.0000001) | (np.abs(data[:,1]) > 0.0000001) ) ]
	
	# Generate plots
	if len(cols) == 2:
		make_2D_density_plot_v1( data[:,0], data[:,1], r'$x$ (fm)', r'$y$ (fm)' )
	else:
		make_2D_density_plot_v2( data[:,0], data[:,1], data[:,2], r'$x$ (fm)', r'$y$ (fm)' )

	# Pause so we can look at the plots	
	#pause()
	

# End of file
