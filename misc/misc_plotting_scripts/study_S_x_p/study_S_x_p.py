#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os, glob

#====================================================
#def pause():
#    programPause = raw_input()

#====================================================
def make_2D_density_plot( xDir, yDir, xLimits, yLimits, \
                            xLabel, yLabel, bws,
                            outputfilename ):
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
                    vmin = 0.0, vmax = 1.0,
                    extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()])

    plt.colorbar(im,fraction=0.026, pad=0.04)

    plt.xlabel(xLabel)
    plt.ylabel(yLabel)

    plt.xlim(xLimits)
    plt.ylim(yLimits)
	
    plt.tight_layout()
    #plt.show()
    fig.savefig(outputfilename)
    print 'Saving to', outputfilename



#====================================================
def generate_plot( data, polarMode, bws, \
                    plotMode, xLimits, yLimits,
                    outputfilename):
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
    elif plotMode == 3:
        xLabel = r'$x_o$ (fm)'
        yLabel = r'$\tau$ (fm)'
    elif plotMode == 4:
        xLabel = r'$x_s$ (fm)'
        yLabel = r'$\tau$ (fm)'
    elif plotMode == 5:
        xLabel = r'$x_o$ (fm)'
        yLabel = r'$t$ (fm$/c$)'
    elif plotMode == 6:
        xLabel = r'$x_s$ (fm)'
        yLabel = r'$t$ (fm$/c$)'
    else:
        print 1/0

    [Kphi,t,x,y,z]=data
	
    # Plot them
    if plotMode == 0:
        if polarMode == 0:
            make_2D_density_plot(z, t, xLimits, yLimits, xLabel, yLabel, bws, outputfilename)
        else:
            safeIndices = np.where(np.greater(t**2, x**2+y**2+z**2))
            tsI = t[safeIndices]
            zsI = z[safeIndices]
            eta = 0.5*np.log( (tsI+zsI)/(tsI-zsI) )
            tau = np.sqrt(tsI**2-zsI**2)
            make_2D_density_plot(eta, tau, xLimits, yLimits, xLabel, yLabel, bws, outputfilename)
    elif plotMode == 1:
        make_2D_density_plot(x, y, xLimits, yLimits, xLabel, yLabel, bws, outputfilename)
    elif plotMode == 2:
        safeIndices = np.where(np.greater(t**2, x**2+y**2+z**2))
        tsI = t[safeIndices]
        xsI = x[safeIndices]
        ysI = y[safeIndices]
        zsI = z[safeIndices]
        r = np.sqrt(xsI**2+ysI**2)
        tau = np.sqrt(tsI**2-zsI**2)
        make_2D_density_plot(r, tau, xLimits, yLimits, xLabel, yLabel, bws, outputfilename)
    elif plotMode >= 3 and plotMode <= 6:
        safeIndices = np.where(np.greater(t**2, x**2+y**2+z**2))
        tsI = t[safeIndices]
        xsI = x[safeIndices]
        ysI = y[safeIndices]
        zsI = z[safeIndices]
        KphisI = Kphi[safeIndices]
        tau = np.sqrt(tsI**2-zsI**2)
        eta = 0.5*np.log((tsI+zsI)/(tsI-zsI))
        if plotMode == 3:
            xo = xsI * np.cos(KphisI) + ysI * np.sin(KphisI)
            make_2D_density_plot(xo, tau, xLimits, yLimits, xLabel, yLabel, bws, outputfilename)
        elif plotMode == 4:
            xs = -xsI * np.sin(KphisI) + ysI * np.cos(KphisI)
            make_2D_density_plot(xs, tau, xLimits, yLimits, xLabel, yLabel, bws, outputfilename)
        elif plotMode == 5:
            xo = xsI * np.cos(KphisI) + ysI * np.sin(KphisI)
            make_2D_density_plot(xo, tsI, xLimits, yLimits, xLabel, yLabel, bws, outputfilename)
        elif plotMode == 6:
            xs = -xsI * np.sin(KphisI) + ysI * np.cos(KphisI)
            make_2D_density_plot(xs, tsI, xLimits, yLimits, xLabel, yLabel, bws, outputfilename)


#====================================================
def generate_plots(filename):
    # Read in name of file from command line
    #filename = sys.argv[1]
    #filename = "C:/Users/Christopher Plumberg/Desktop/Research/Lund"\
    #            +"/Multiplicity_dependence_of_HBT_w_Pythia"\
    #            +"/study_S_x_p_results_KLmax0.01/S_x_p_N1_11_0.0_0.1.dat"

    # Load file
    data = np.loadtxt(filename).T
    
    KphiStem=''
    #xSector=(sys.argv[2]=="True")
    #ySector=(sys.argv[3]=="True")
    xSector=False
    ySector=False
    if xSector:
        data=data.T
        data[np.where(data[:,0] > np.pi),0] -= 2.0*np.pi
        data = data[np.where(np.abs(data[:,0])<0.125*np.pi)]
        data=data.T
        KphiStem = '_xPos'
    elif ySector:
        data=data.T
        data = data[np.where(np.abs(data[:,0]-0.5*np.pi)<0.125*np.pi)]
        data=data.T
        KphiStem = '_yPos'
    
    filenameStem = (os.path.splitext(filename)[0]).replace(".","")
	
    # Generate plots
    generate_plot( data, 0, [0.1, 0.1], 0, [-10.0, 10.0], [0.0, 15.0], filenameStem+KphiStem+'_z_t.pdf' )
    generate_plot( data, 1, [0.1, 0.1], 0, [-5.0, 5.0], [0.0, 5.0], filenameStem+KphiStem+'_eta_tau.pdf' )
    generate_plot( data, 0, [0.1, 0.1], 1, [-2.5, 2.5], [-2.5, 2.5], filenameStem+KphiStem+'_x_y.pdf' )
    generate_plot( data, None, [0.1, 0.1], 2, [0.0, 5.0], [0.0, 5.0], filenameStem+KphiStem+'_r_tau.pdf' )
    generate_plot( data, None, [0.1, 0.1], 3, [-5.0, 5.0], [0.0, 5.0], filenameStem+KphiStem+'_xo_tau.pdf' )
    generate_plot( data, None, [0.1, 0.1], 4, [-5.0, 5.0], [0.0, 5.0], filenameStem+KphiStem+'_xs_tau.pdf' )
    generate_plot( data, None, [0.1, 0.1], 5, [-5.0, 5.0], [0.0, 5.0], filenameStem+KphiStem+'_xo_t.pdf' )
    generate_plot( data, None, [0.1, 0.1], 6, [-5.0, 5.0], [0.0, 5.0], filenameStem+KphiStem+'_xs_t.pdf' )



#====================================================
if __name__ == "__main__":
    
    for filename in glob.glob("C:/Users/Christopher Plumberg/Desktop/Research/Lund" \
                                + "/Multiplicity_dependence_of_HBT_w_Pythia" \
                                + "/study_S_x_p_results_KLmax0_01/S_x_p_N*.dat"):
        generate_plots(filename)

# End of file
