#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

# Useful constants, settings, etc.
mmTOfm          = 1.e12
'''
includeThermal  = True
includeDecays   = True
usingLogDensity = True
numberOfBins    = 1000

# 0 means plot (z,t), 1 means plot (x,y)
plotMode        = 1

# Plotting settings ( x and y refer to figure axes,
#					  not necessarily physical coordinates )
xLimits = [-15.0,15.0]
yLimits = [-15.0,15.0]
'''

#====================================================
def pause():
    programPause = raw_input()


#====================================================
def load_OSCAR_file(filename):
	#print 'Reading in '+filename+'...'
	with open(filename,"r") as f:
		all_data=[x.split() for x in f.readlines()]

	#print 'Getting event records...'
	eventRecords=[[int(x[0]), all_data.index(x), int(x[1])] for x in all_data if len(x)==2]
	#all_data=[]

	#print 'Splitting event information...'
	split_data=[ all_data[ ER[1]+1:ER[1]+ER[2]+1 ] for ER in eventRecords ]
	#eventRecords=[]

	#print 'Processing...'
	ensemble=[[map( (lambda x: float(x)) , col2) for col2 in col1] for col1 in split_data]
	#split_data=[]

	#print 'Check: len(ensemble) =', len(ensemble)
	#print 'Check: len(ensemble[0]) =', len(ensemble[0])
	#print ensemble[0][0:10]

	#data = np.array([particle for event in ensemble[0:1] for particle in event])
	data = np.array([particle for event in ensemble for particle in event])
	#ensemble=[]

	#print 'Finished! data.shape =', data.shape
	return data




#====================================================
def make_2D_density_plot( xDir, yDir, xLimits, yLimits, \
							xLabel, yLabel, nbins, logDensity ):
	fig, ax = plt.subplots()

	H, xedges, yedges = np.histogram2d(xDir, yDir, bins=nbins)

	xcenters = (xedges[:-1] + xedges[1:]) / 2
	ycenters = (yedges[:-1] + yedges[1:]) / 2

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
def generate_plot( data, includeThermal, includeDecays, \
					usingLogDensity, numberOfBins, \
					plotMode, xLimits, yLimits ):
	# Set the quantities to plot
	[t, x, y, z] = mmTOfm * data[:,[8,9,10,11]].T
	resID = data[:,1]
	
	# Take relevant slices of selected quantities
	[xLabel, yLabel] = [None, None]
	[xCondition, yCondition, zCondition, tCondition] \
		= [True, True, True, True]
	if plotMode == 0:
		xLabel = r'$z$ (fm)'
		yLabel = r'$t$ (fm/$c$)'
		zCondition = (z > xLimits[0]) & (z < xLimits[1])
		tCondition = (t > yLimits[0]) & (t < yLimits[1])
	elif plotMode == 1:
		xLabel = r'$x$ (fm)'
		yLabel = r'$y$ (fm)'
		xCondition = (x > xLimits[0]) & (x < xLimits[1])
		yCondition = (y > yLimits[0]) & (y < yLimits[1])
	else:
		print 1/0
	thermalCondition = includeThermal & (resID==0)
	decayCondition = includeDecays & (resID==1)
	chosenSlice = np.where( tCondition & xCondition & yCondition & zCondition & ( thermalCondition | decayCondition ) )
	[t,x,y,z]=np.c_[t,x,y,z][chosenSlice].T
	
	# Plot them
	if plotMode == 0:
		make_2D_density_plot(z, t, xLimits, yLimits, xLabel, yLabel, numberOfBins, usingLogDensity)
	elif plotMode == 1:
		make_2D_density_plot(x, y, xLimits, yLimits, xLabel, yLabel, numberOfBins, usingLogDensity)



#====================================================
if __name__ == "__main__":
	# Read in name of file from command line
	filename = sys.argv[1]
	
	# Load file in OSCAR-type format
	data = load_OSCAR_file(filename)
	
	# Generate plots
	#generate_plot( data, True, True, True, 50, 0, [-15.0, 15.0], [-1.0, 20.0] )
	#generate_plot( data, True, False, True, 1000, 0, [-15.0, 15.0], [-1.0, 20.0] )
	#generate_plot( data, False, True, True, 1000, 0, [-15.0, 15.0], [-1.0, 20.0] )
	generate_plot( data, True, True, False, 50, 1, [-5.0, 5.0], [-5.0, 5.0] )
	#generate_plot( data, True, False, False, 100, 1, [-5.0, 5.0], [-5.0, 5.0] )
	#generate_plot( data, False, True, False, 100, 1, [-5.0, 5.0], [-5.0, 5.0] )

	# Pause so we can look at the plots	
	#pause()
	

# End of file
