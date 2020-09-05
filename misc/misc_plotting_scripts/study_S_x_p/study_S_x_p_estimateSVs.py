#!/usr/bin/env python
#====================

#from scipy import stats
import scipy.special as sp, scipy.stats
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
    hist_dist = stats.rv_histogram((sums, bins))
    #pars = gamma.fit(weightedData, floc=0.0)
    pars = dist.fit(hist_dist.rvs(size=10000000), floc=0.0)
    a1, loc1, scale1 = pars
    #print minimum, maximum, a1, loc1, scale1
    return a1, scale1

def estimate_SVs(data):
    # parse data and impose limits (in fm)
    print 'Loading and parsing data...'
    xmin, xmax = -10, 10
    ymin, ymax = -10, 10
    
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
    print '-------------------------------------'
    print 'Computing variances and radii...'
    print '-------------------------------------'
    print '  --> means:'
    mean_xo, mean_xs, mean_xl, mean_t = np.mean(xo), np.mean(xs), np.mean(xl), np.mean(t)
    print '  --> ', mean_xo, mean_xs, mean_xl, mean_t
    var_xo, var_xs, var_xl, var_t = np.var(xo), np.var(xs), np.var(xl), np.var(t)
    cov_xo_t = np.mean(xo*t)-np.mean(xo)*np.mean(t)
    print '-------------------------------------'
    print '  --> (co)variances:'
    print '  --> ', var_xo, var_xs, var_xl, var_t
    print '  --> ', cov_xo_t, np.cov(xo, t, bias=True)[0][1], np.cov(xo, t, bias=False)[0][1]
    print '-------------------------------------'
    print '  --> radii:'
    R2o, R2s, R2l = var_xo - 2.0*betaT*(np.mean(xo*t)-np.mean(xo)*np.mean(t))+betaT**2*np.var(t),\
          np.var(xs), np.var(xl)
    print '  --> ', R2o, R2s, R2l 
          
    # Estimate t and z variances by imposing cut on tau
    tauMax = 10.0
    tau, eta = np.sqrt(t**2-z**2), 0.5*np.log((t+z)/(t-z))
    milne = np.c_[ tau, eta ]
    milne = milne[np.where(milne[:,0]<tauMax)]
    tTrunc = milne[:,0]*np.cosh(milne[:,1])
    zTrunc = milne[:,0]*np.sinh(milne[:,1])
    print '-------------------------------------'
    print '  --> truncated t/z means/variances:'
    print '  --> ', np.mean(zTrunc), np.mean(tTrunc), np.var(zTrunc), np.var(tTrunc)
    
    # Estimate t and z variances instead from fitting to Gamma distribution
    print '-------------------------------------'
    print '  --> MLE t/z means/variances:'
    tau, eta = np.sqrt(t**2-z**2), 0.5*np.log((t+z)/(t-z))
    # N.B. - extra factor of tau pulled out!
    alpha, beta = do_MLE_withWeights(tau, scipy.stats.gamma, 0.01)
    tau3mom = scipy.stats.gamma.expect(lambda TAU: TAU**3, args=(alpha,), loc=0, scale=beta)
    tau2mom = scipy.stats.gamma.expect(lambda TAU: TAU**2, args=(alpha,), loc=0, scale=beta)
    tau1mom = scipy.stats.gamma.expect(lambda TAU: TAU, args=(alpha,), loc=0, scale=beta)
    cosh_eta_mom = np.mean(np.cosh(eta))
    sinh_eta_mom = np.mean(np.sinh(eta))
    cosh2_eta_mom = np.mean(np.cosh(eta)**2)
    sinh2_eta_mom = np.mean(np.sinh(eta)**2)
    t2_expect = tau3mom * cosh2_eta_mom / tau1mom
    z2_expect = tau3mom * sinh2_eta_mom / tau1mom
    t_expect = tau2mom * cosh_eta_mom / tau1mom
    z_expect = tau2mom * sinh_eta_mom / tau1mom
    print '  --> ', z_expect, t_expect, z2_expect - z_expect**2, t2_expect - t_expect**2
        
                
#====================================================
if __name__ == "__main__":
    # Read in name of file from command line
    filename = sys.argv[1]
    KT = float(sys.argv[2])
	
    # Load file
    data = np.loadtxt(filename).T
    
    estimate_SVs(data)
    
    # Scratch and backup below this line
    '''xo=xo[np.argsort(np.abs(xo))]
    
    #print np.min(xo), np.max(xo)
    
    for i in range(len(xo),2,-1):
        print float(i)/float(len(xo)), np.var(xo[:i]), np.min(xo[:i]), np.max(xo[:i])'''
    
    '''print '-------------------------------------'
    print 'Regular stats:'
    print '<x_o>, <x_s>, <x_l>, <t>:'
    print np.mean(xo), np.mean(xs), np.mean(xl), np.mean(t)
    print
    print '<(x_o-<x_o>)^2>, <(x_s-<x_s>)^2>, <(x_l-<x_l>)^2>, <(x_o-<x_o>)(t-<t>)>, <(t-<t>)^2>:'
    print np.var(xo), np.var(xs), np.var(xl), \
            np.mean(xo*t)-np.mean(xo)*np.mean(t), np.var(t)
    print
    print 'R^2_o, R^2_s, R^2_l:'
    print np.var(xo) - 2.0*betaT*(np.mean(xo*t)-np.mean(xo)*np.mean(t))+betaT**2*np.var(t),\
          np.var(xs), np.var(xl)

    print '-------------------------------------'
    print 'Trimmed stats:'
    for ul in [10, 100, 1000, 10000, 100000]:
        print 'ul =', ul
        print '<x_o>, <x_s>, <x_l>, <t>:'
        print stats.tmean(xo, (-ul, ul)), stats.tmean(xs, (-ul, ul)), \
               stats.tmean(xl, (-ul, ul)), stats.tmean(t, (None, ul))
        print '<(x_o-<x_o>)^2>, <(x_s-<x_s>)^2>, <(x_l-<x_l>)^2>, <(t-<t>)^2>:'
        print stats.tvar(xo, (-ul, ul)), stats.tvar(xs, (-ul, ul)), \
               stats.tvar(xl, (-ul, ul)), stats.tvar(t, (None, ul))'''

# End of file