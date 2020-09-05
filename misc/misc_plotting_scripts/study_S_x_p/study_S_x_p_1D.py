#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from scipy.stats import gamma, nakagami, rv_histogram
from scipy import optimize
        
def GammaDistribution(pars, x):
    #alpha = 8.03748
    #beta = 0.18678
    [alpha, beta] = pars
    return np.exp(-x/beta)*x**(alpha-1.0)/(beta**alpha * sp.gamma(alpha))
    
def NakagamiDistribution(pars, x):
    [m, w] = pars
    return 2.0*(m/w)**m * x**(2.0*m-1.0) * np.exp(-m*x**2/w)/sp.gamma(m)


bounds = [(0, 10), (0, 10), (0, 10)]
def optimizer(func, x0, args, disp):
    res = optimize.differential_evolution(func, bounds, args, tol=1e-1)
    return res.x

def do_MLE(data, minimum, maximum):
    pars = gamma.fit(data[np.where((minimum<=data)&(maximum>=data))], floc=0.0)
    a1, loc1, scale1 = pars
    #print minimum, maximum, a1, loc1, scale1, \
    #      gamma.nnlf(pars, data[np.where((minimum<=data)&(maximum>=data))]), \
    #      a1 * scale1, a1 * scale1**2
    return a1, scale1


def do_MLE_withWeights(data, dist, minimum, maximum, bw):
    nbins = int( (maximum - minimum) / bw )
    rawData = data[np.where((minimum<=data)&(maximum>=data))]
    sums, bins = np.histogram(rawData, bins=nbins, range=[minimum, maximum])
    bincenters = (lambda v: 0.5*(v[1:]+v[:-1]))(bins)
    sums, bins = np.histogram(bins[:-1], bins=nbins, range=[minimum, maximum], \
                              density=True, weights=sums/bincenters)
    hist_dist = rv_histogram((sums, bins))
    #pars = gamma.fit(weightedData, floc=0.0)
    pars = dist.fit(hist_dist.rvs(size=10000000), floc=0.0)
    a1, loc1, scale1 = pars
    print minimum, maximum, a1, loc1, scale1, a1 * scale1, a1 * scale1**2
    return a1, scale1


#====================================================
def make_1D_density_plot( data, limits, bw ):
    nbins = int( (limits[1] - limits[0]) / bw )
    fig, ax = plt.subplots()

    #bincenters = (lambda v: 0.5*(v[1:]+v[:-1]))(np.arange(limits[0], limits[1]+bw, bw))
    #print len(bincenters), nbins, data.shape
    n, bins, patches = plt.hist( data, nbins, range=limits, density=True, \
                                 facecolor='g', alpha=0.75 )
    #print np.sum(0.1*n)
    #print 1/0
    #counts, bins = np.histogram(data, bins=nbins, range=limits, density=True )
    #print len(bins), len(counts)
    #print bins
    #bincenters = (lambda v: 0.5*(v[1:]+v[:-1]))(bins)
    #plt.hist(bins[:-1], bins, density=True, weights=counts/bincenters, \
    #         facecolor='g', alpha=0.75)
    #plt.hist(bins[:-1], bins, density=True, weights=counts, \
    #         facecolor='g', alpha=0.75)
        
    '''eps=0
    x = np.linspace(0.0, 4.0, 1001)
    #plt.plot(x, GammaDistribution(do_MLE(data, eps, 3.0), x), 'r-')
    #plt.plot(x, GammaDistribution(do_MLE(data, eps, 4.0), x), 'b-')
    #plt.plot(x, GammaDistribution(do_MLE(data, eps, 5.0), x), '-', color='purple')
    #plt.plot(x, NakagamiDistribution([2.03196, 2.58968], x), '-', color='orange')
    plt.plot(x, GammaDistribution(do_MLE_withWeights(data, gamma, eps, 3.0, 0.01), x), 'r-')
    plt.plot(x, GammaDistribution(do_MLE_withWeights(data, gamma, eps, 4.0, 0.01), x), 'b-')
    plt.plot(x, GammaDistribution(do_MLE_withWeights(data, gamma, eps, 5.0, 0.01), x), '-', color='purple')
    plt.plot(x, GammaDistribution(do_MLE_withWeights(data, gamma, eps, 100000.0, 0.01), x), '-', color='magenta')
    plt.plot(x, NakagamiDistribution([1.76517, 1.96998], x), '-', color='orange')
    plt.plot(x, NakagamiDistribution([2.03196, 2.58968], x), '-', color='yellow')
    plt.plot(x, NakagamiDistribution([1.50433, 2.10114], x), '-', color='cyan')'''

    x = np.linspace(limits[0], limits[1], 1001)
    s = 1.15954
    y = np.exp(-0.5*x**2/s**2)/(np.sqrt(2.0*np.pi)*s)
    plt.plot(x, y, 'r-')
    
    #plt.xlim(limits)
    plt.tight_layout()
    plt.show()



#====================================================
if __name__ == "__main__":
    # Read in name of file from command line
    #filename = sys.argv[1]
    filename = "C:/Users/Christopher Plumberg/Desktop/Research/Lund"\
                +"/Multiplicity_dependence_of_HBT_w_Pythia/Figures"\
                +"/study_S_x_p_results_KLmax0.01/S_x_p_N1_11_0.0_0.1.dat"
	
    # Load file
    data = np.loadtxt(filename).T
    
    [Kphi,t,x,y,z]=data
    safeIndices = np.where(np.greater(t**2, x**2+y**2+z**2))
    tsI = t[safeIndices]
    zsI = z[safeIndices]
    eta = 0.5*np.log( (tsI+zsI)/(tsI-zsI) )
    tau = np.sqrt(tsI**2-zsI**2)
    make_1D_density_plot( eta, [-8.0, 8.0], 0.1 )
    make_1D_density_plot( tau, [0.0, 4.0], 0.1 )
    
#End of file