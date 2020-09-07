#!/usr/bin/env python
#====================

import scipy.special as sp
from scipy.stats import gamma, rv_histogram
import numpy as np
import sys

# PDF of Gamma distribution
def GammaDistribution(pars, x):
    [alpha, beta] = pars
    return np.exp(-x/beta)*x**(alpha-1.0)/(beta**alpha * sp.gamma(alpha))

# PDF of Normal distribution
def NormalDistribution(pars, x):
    [mu, sigma] = pars
    return np.exp(-0.5*(x-mu)**2/sigma**2)/(sigma*np.sqrt(2.0*np.pi))

# Maximum likelihood estimation of distribution parameters
def do_MLE_withWeights(data, dist, bw):
    minimum, maximum = np.min(data), np.max(data)
    nbins = int( (maximum - minimum) / bw )
    rawData = data
    sums, bins = np.histogram(rawData, bins=nbins, range=[minimum, maximum])
    bincenters = (lambda v: 0.5*(v[1:]+v[:-1]))(bins)
    sums, bins = np.histogram(bins[:-1], bins=nbins, range=[minimum, maximum], \
                              density=True, weights=sums/bincenters)
    hist_dist = rv_histogram((sums, bins))
    #pars = gamma.fit(weightedData, floc=0.0)
    pars = dist.fit(hist_dist.rvs(size=10000000), floc=0.0)
    a1, loc1, scale1 = pars
    #print minimum, maximum, a1, loc1, scale1
    return a1, scale1

def estimate_SVs(data, Nmin, Nmax, KTmin, KTmax, verbose=False):
    # parse data and impose limits (in fm)
    if verbose:
        print 'Loading and parsing data...'
    xmin, xmax = -10, 10
    ymin, ymax = -10, 10
    
    KT = 0.5*(float(KTmin)+float(KTmax))
    
    data = data.T
    xCondition = (xmin < data[:,2]) & (xmax > data[:,2])
    yCondition = (ymin < data[:,3]) & (ymax > data[:,3])
    causalCondition = (data[:,1]**2 > data[:,2]**2+data[:,3]**2+data[:,4]**2)
    data = data[np.where((xCondition) & (yCondition) & (causalCondition))].T
    [Kphi,t,x,y,z]=data
    
    xo = x*np.cos(Kphi) + y*np.sin(Kphi)
    xs = -x*np.sin(Kphi) + y*np.cos(Kphi)
    xl = z
    
    mpi = 0.13957
    betaT = KT / np.sqrt(mpi**2 + KT**2)

    # compute moments, etc. directly
    mean_xo, mean_xs, mean_xl, mean_t = np.mean(xo), np.mean(xs), np.mean(xl), np.mean(t)
    var_xo, var_xs, var_xl, var_t = np.var(xo), np.var(xs), np.var(xl), np.var(t)
    cov_xo_t = np.mean(xo*t)-np.mean(xo)*np.mean(t)
    R2o, R2s, R2l = var_xo \
                    - 2.0*betaT*(np.mean(xo*t) - np.mean(xo)*np.mean(t)) \
                    + betaT**2*np.var(t),\
                    np.var(xs), np.var(xl)
    if verbose:
        print '-------------------------------------'
        print 'Computing variances and radii...'
        print '-------------------------------------'
        print '  --> means:'
        print '  --> ', mean_xo, mean_xs, mean_xl, mean_t
        print '-------------------------------------'
        print '  --> (co)variances:'
        print '  --> ', var_xo, var_xs, var_xl, var_t
        print '  --> ', cov_xo_t, np.cov(xo, t, bias=True)[0][1], np.cov(xo, t, bias=False)[0][1]
        print '-------------------------------------'
        print '  --> radii:'
        print '  --> ', R2o, R2s, R2l 
          
    # Estimate t and z variances by imposing cut on tau
    tauMax = 10.0
    tau, eta = np.sqrt(t**2-z**2), 0.5*np.log((t+z)/(t-z))
    milne = np.c_[ tau, eta, xo, xs, z, t ]
    milne = milne[np.where(milne[:,0]<tauMax)]
    tTrunc = milne[:,0]*np.cosh(milne[:,1])
    zTrunc = milne[:,0]*np.sinh(milne[:,1])
    truncR2o, truncR2s, truncR2l = var_xo \
                                   - 2.0*betaT*np.cov(milne[:,2], milne[:,5], bias=True)[0][1]\
                                   + betaT**2*np.var(tTrunc),\
                                   np.var(xs), np.var(zTrunc)
    if verbose:
        print '-------------------------------------'
        print '  --> truncated t/z means/variances:'
        print '  --> ', np.mean(zTrunc), np.mean(tTrunc), np.var(zTrunc), np.var(tTrunc)
        print '  --> ', np.mean(milne[:,4]), np.mean(milne[:,5]), np.var(milne[:,4]), np.var(milne[:,5])
        print '  --> ', np.cov(milne[:,2], milne[:,5], bias=True)[0][1]
        print '-------------------------------------'
        print '  --> truncated radii:'
        print '  --> ', truncR2o, truncR2s, truncR2l
        
    return var_xo, np.cov(milne[:,2], milne[:,5], bias=True)[0][1], np.var(tTrunc),\
           np.var(xs), np.var(zTrunc), truncR2o, truncR2s, truncR2l
    

                
#====================================================
if __name__ == "__main__":
    # Read in name of file from command line
    # multiplicites and KT values
    multmins = ['1','12','17','23','29','35','42','52']
    multmaxs = ['11','16','22','28','34','41','51','151']
    KTmins = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7']
    KTmaxs = ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8']
    
    allResults = []
    for iMult, Nmin in enumerate(multmins):
        for iKT, KTmin in enumerate(KTmins):
            Nmax, KTmax = multmaxs[iMult], KTmaxs[iKT]
            path = sys.argv[1]
            filename = path + '/' + 'S_x_p_N' + Nmin + '_' + Nmax + '_' + KTmin + '_' + KTmax + '.dat'
            print 'Processing', filename
        
            # Load file
            data = np.loadtxt(filename).T
            
            SvsAndTruncRadii = estimate_SVs(data, Nmin, Nmax, KTmin, KTmax)
            resultsForThisFile = [Nmin, Nmax, KTmin, KTmax] + list(SvsAndTruncRadii)
            allResults.append(resultsForThisFile)
            outfilename = path + '/' + 'SVsAndRadii_N' + Nmin + '_' + Nmax + '_' + KTmin + '_' + KTmax + '.dat'
            np.savetxt(outfilename, np.array(resultsForThisFile), fmt='%s')
            

    outfilename = path + '/' + 'SVsAndRadii_all.dat'
    np.savetxt(outfilename, np.array(allResults), fmt='%s')

# End of file