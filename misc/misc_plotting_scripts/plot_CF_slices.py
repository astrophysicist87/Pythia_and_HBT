#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate
#from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.ndimage import gaussian_filter
import sys, glob

#====================================================
GeVtoMeV = 1000.0
showInsteadOfSave = True

#filenames = sys.argv[1:]
filenames = ['thermal_N10000000_enhanceMode0_fixedQRef_CF.dat', \
				#'thermal_N1000000_enhanceMode0_fixedQRef_CF.dat', \
				'thermal_N1000000_enhanceMode1_CF.dat', \
				'thermal_N1000000_enhanceMode0_noBE_CF.dat']

labels = [r'$1+\exp\left( -Q^2\!R^2\!\right)$, $N_{\mathrm{ev}}=10^7$', \
			#r'$1+\exp\left( -Q^2 R^2 \right)$, $N_{\mathrm{ev}}=10^6$', \
			r'$1+j_0\left( Q \left| \Delta \vec{x} \right| \right)$, $N_{\mathrm{ev}}=10^6$', \
			r'Vertex tracking, $N_{\mathrm{ev}}=10^6$']

colors = ['blue', 'red', 'green', 'purple', 'orange', 'yellow']

#====================================================
def pause():
    programPause = raw_input("Press ENTER to continue...")
#====================================================
   
'''def stack(arrays, axis=0):
	arrays = [np.asanyarray(arr) for arr in arrays]
	if not arrays:
		raise ValueError('need at least one array to stack')
	
	shapes = set(arr.shape for arr in arrays)
	if len(shapes) != 1:
		raise ValueError('all input arrays must have the same shape')
	
	result_ndim = arrays[0].ndim + 1
	if not -result_ndim <= axis < result_ndim:
		msg = 'axis {0} out of bounds [-{1}, {1})'.format(axis, result_ndim)
		raise IndexError(msg)
	if axis < 0:
		axis += result_ndim
	
	sl = (slice(None),) * axis + (_nx.newaxis,)
	expanded_arrays = [arr[sl] for arr in arrays]
	return _nx.concatenate(expanded_arrays, axis=axis)'''

#====================================================
#====================================================
def plot_out_slice(data, KT):
	fig, ax = plt.subplots(1, 1, figsize=(15,10))
	plotfontsize = 20

	i=0
	for dataThisFile in data:
		thisData = dataThisFile[np.where((dataThisFile[:,0] == KT) & (dataThisFile[:,4] == 0) & (dataThisFile[:,5] == 0))]
		ax.plot(thisData[:,3], thisData[:,10], color=colors[i], label=labels[i])
		ax.plot(thisData[:,3], 1.0 + 0.0*thisData[:,10], color='black')
		
		# include error bands
		plt.fill_between(thisData[:,3],\
							thisData[:,10]-thisData[:,11],\
							thisData[:,10]+thisData[:,11],\
							alpha=0.25, edgecolor=colors[i],\
							facecolor=colors[i])

		xmin, xmax = -0.25, 0.25
		ymin, ymax = 0.95, 2.025
		lims = [xmin, xmax, ymin, ymax]
		plt.axis(lims)
		
		ax.tick_params(axis='both', labelsize=plotfontsize)
		ax.set_xlabel(r'$q_o$ (GeV)', fontsize=plotfontsize+8)
		ax.set_ylabel(r'$C(q_o)$', fontsize=plotfontsize+8)
		#ax.legend(loc='best', fontsize=plotfontsize)

		i+=1


	if showInsteadOfSave:
		plt.show(block=False)
	else:
		outfilename = './CFs_vs_qo.pdf'
		print 'Saving to', outfilename
		plt.savefig(outfilename, bbox_inches='tight', dpi=250)
		plt.close()

#====================================================
#====================================================
def plot_side_slice(data, KT):
	fig, ax = plt.subplots(1, 1, figsize=(15,10))
	plotfontsize = 20

	i=0
	for dataThisFile in data:
		thisData = dataThisFile[np.where((dataThisFile[:,0] == KT) & (dataThisFile[:,3] == 0) & (dataThisFile[:,5] == 0))]
		ax.plot(thisData[:,4], thisData[:,10], color=colors[i], label=labels[i])
		ax.plot(thisData[:,4], 1.0 + 0.0*thisData[:,10], color='black')
		
		# include error bands
		plt.fill_between(thisData[:,4],\
							thisData[:,10]-thisData[:,11],\
							thisData[:,10]+thisData[:,11],\
							alpha=0.25, edgecolor=colors[i],\
							facecolor=colors[i])

		xmin, xmax = -0.25, 0.25
		ymin, ymax = 0.95, 2.025
		lims = [xmin, xmax, ymin, ymax]
		plt.axis(lims)
		
		ax.tick_params(axis='both', labelsize=plotfontsize)
		ax.set_xlabel(r'$q_s$ (GeV)', fontsize=plotfontsize+8)
		ax.set_ylabel(r'$C(q_s)$', fontsize=plotfontsize+8)
		#ax.legend(loc='best', fontsize=plotfontsize)

		i+=1


	if showInsteadOfSave:
		plt.show(block=False)
	else:
		outfilename = './CFs_vs_qs.pdf'
		print 'Saving to', outfilename
		plt.savefig(outfilename, bbox_inches='tight', dpi=250)
		plt.close()

#====================================================
#====================================================
def plot_long_slice(data, KT):
	fig, ax = plt.subplots(1, 1, figsize=(15,10))
	plotfontsize = 20

	i=0
	for dataThisFile in data:
		thisData = dataThisFile[np.where((dataThisFile[:,0] == KT) & (dataThisFile[:,3] == 0) & (dataThisFile[:,4] == 0))]
		ax.plot(thisData[:,5], thisData[:,10], color=colors[i], label=labels[i])
		ax.plot(thisData[:,5], 1.0 + 0.0*thisData[:,10], color='black')
		
		# include error bands
		plt.fill_between(thisData[:,5],\
							thisData[:,10]-thisData[:,11],\
							thisData[:,10]+thisData[:,11],\
							alpha=0.25, edgecolor=colors[i],\
							facecolor=colors[i])

		xmin, xmax = -0.25, 0.25
		ymin, ymax = 0.95, 2.025
		lims = [xmin, xmax, ymin, ymax]
		plt.axis(lims)
		
		ax.tick_params(axis='both', labelsize=plotfontsize)
		ax.set_xlabel(r'$q_l$ (GeV)', fontsize=plotfontsize+8)
		ax.set_ylabel(r'$C(q_l)$', fontsize=plotfontsize+8)
		ax.legend(loc='best', fontsize=plotfontsize+2)

		i+=1


	if showInsteadOfSave:
		plt.show(block=False)
	else:
		outfilename = './CFs_vs_ql.pdf'
		print 'Saving to', outfilename
		plt.savefig(outfilename, bbox_inches='tight', dpi=250)
		plt.close()


#====================================================
if __name__ == "__main__":
	# Load files into single data object
	#print filenames
	data = np.array([np.loadtxt(f) for f in filenames])

	#print data.shape
	plot_out_slice(data, 0.1)
	plot_side_slice(data, 0.1)
	plot_long_slice(data, 0.1)

	if showInsteadOfSave:
		pause()



# End of file
