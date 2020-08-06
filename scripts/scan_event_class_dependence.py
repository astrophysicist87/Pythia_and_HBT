import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys, os

GeVToMeV = 1000.0
recycle = True

columnLabels = ['dNdeta','KT',\
                'lambda','R2o','R2s','R2l',\
                'err(lambda)','err(R2o)','err(R2s)','err(R2l)']
#multLabels = [r'$N_{ch}=1-11$', r'$N_{ch}=12-16$', r'$N_{ch}=17-22$', \
#              r'$N_{ch}=23-28$', r'$N_{ch}=29-34$', r'$N_{ch}=35-41$', \
#              r'$N_{ch}=42-51$', r'$N_{ch}=52-151$', r'$N_{ch}=152-\infty$']
multLabels = [r'$0-100\%$']
#KTLabels = [r'$K_T=0-100$ MeV', r'$K_T=100-200$ MeV', r'$K_T=200-300$ MeV', \
#            r'$K_T=300-400$ MeV', r'$K_T=400-500$ MeV', r'$K_T=500-600$ MeV', \
#            r'$K_T=600-700$ MeV', r'$K_T=700-800$ MeV']
KTLabels = [r'$K_T=0-100$ MeV', r'$K_T=100-200$ MeV', r'$K_T=200-300$ MeV', \
            r'$K_T=300-400$ MeV', r'$K_T=400-500$ MeV', r'$K_T=500-600$ MeV']

R2iLabels = dict({'R2o':r'$R_o$ $($fm$)$','R2s':r'$R_s$ $($fm$)$','R2l':r'$R_l$ $($fm$)$'})

cols = dict(zip(columnLabels,range(len(columnLabels))))
CFcols=dict({'out': 0,'side': 1,'long': 2})
colors = ['blue','red','green','purple','cyan','magenta','orange','yellow','teal']
styles = ['o','*','v','^','s','D','P','X','d']

def pause():
    programPause = raw_input("")

#def get_event_class_mean_dNchdeta(filename):
#    return np.mean(np.loadtxt(filename, usecols=4))
def get_event_class_mean_dNchdeta(directory):
    if os.path.exists(directory + '/CF_results/mean_dNchdeta.dat') and recycle:
        return np.loadtxt(directory + '/CF_results/mean_dNchdeta.dat')
    else:
        result = np.mean(np.loadtxt(directory \
                                    + '/CF_results/event_class_multiplicities.dat',\
                                    usecols=6))
        np.savetxt(directory + '/CF_results/mean_dNchdeta.dat', np.array([result]))
        return result

def get_dNchdeta_and_R2i_data(directory):
    mean_dNchdeta = get_event_class_mean_dNchdeta(directory)
    #mean_dNchdeta = get_event_class_mean_dNchdeta(directory + '/CF_results/event_class_multiplicities.dat')
    R2i_vs_KT = np.loadtxt(directory + '/fit_results/HBTradii.dat', \
                           usecols=(0,3,4,5,6,10,11,12,13))[0:6]
    return np.c_[ np.broadcast_to( mean_dNchdeta, len(R2i_vs_KT) ), R2i_vs_KT ]
    
def get_dNchdeta_and_CF_data(directory):
    nMult = len(multLabels)
    #nq = 31
    nq = 51
    #nq = 63
    nKT = len(KTLabels)
    mean_dNchdeta = get_event_class_mean_dNchdeta(directory)
    CF_vs_KT = np.loadtxt(directory + '/CF_results/HBT_pipiCF.dat', usecols=(0,3,4,5,10,11))
    CF_vs_KT_o = CF_vs_KT[np.where((CF_vs_KT[:,2]==0) & (CF_vs_KT[:,3]==0))].reshape([8,nq,nKT])
    #print CF_vs_KT_o[:,:,[0,1,4,5]].shape
    CF_vs_KT_o = np.c_[ np.broadcast_to( mean_dNchdeta, [8,nq,1] ), CF_vs_KT_o[:,:,[0,1,4,5]] ]
    CF_vs_KT_s = CF_vs_KT[np.where((CF_vs_KT[:,1]==0) & (CF_vs_KT[:,3]==0))].reshape([8,nq,nKT])
    CF_vs_KT_s = np.c_[ np.broadcast_to( mean_dNchdeta, [8,nq,1] ), CF_vs_KT_s[:,:,[0,2,4,5]] ]
    CF_vs_KT_l = CF_vs_KT[np.where((CF_vs_KT[:,1]==0) & (CF_vs_KT[:,2]==0))].reshape([8,nq,nKT])
    CF_vs_KT_l = np.c_[ np.broadcast_to( mean_dNchdeta, [8,nq,1] ), CF_vs_KT_l[:,:,[0,3,4,5]] ]
    return np.stack( (CF_vs_KT_o, CF_vs_KT_s, CF_vs_KT_l) )
    


def plot_R2i_vs_KT(data, R2i):
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linestyle='-', linewidth=1.0)
    for multIndex in range(len(data)):
        dataSlice=data[multIndex]
        R2iVec=dataSlice[:,cols[R2i]]
        R2iVec[R2iVec<0] = 0.0
        ax.plot(dataSlice[:,cols['KT']]*GeVToMeV, \
                np.sqrt(R2iVec), \
                '-'+styles[multIndex], \
                color=colors[multIndex])
    ax.set_xlabel(r'$K_T$ (MeV)', fontsize=16)
    ax.set_ylabel(R2iLabels[R2i], fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    plt.tight_layout()
    #plt.show(block=False)
    #plt.show()
    fig.savefig("./"+R2i+"_vs_KT.pdf")
    return None

def plot_R2i_vs_KT_dummy(data):
    fig, ax = plt.subplots()
    for multIndex in range(len(data)):
        dataSlice=data[multIndex]
        ax.plot(dataSlice[:,cols['KT']]*GeVToMeV, \
                -1.0+0.0*dataSlice[:,cols['KT']], \
                '-'+styles[multIndex], \
                color=colors[multIndex], \
                label=multLabels[multIndex])
    ax.set_ylim(bottom=0.0)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.legend(ncol=1, loc='center', fontsize=18)
    plt.tight_layout()
    fig.savefig("./R2i_vs_KT_dummy.pdf")
    #plt.show(block=False)
    #plt.show()
    return None

def plot_R2i_vs_mult(data, R2i):
    fig, ax = plt.subplots()
    ax.axhline(0, color='k', linestyle='-', linewidth=1.0)
    for KTIndex in range(data.shape[1]):
        dataSlice=data[:,KTIndex,:]
        R2iVec=dataSlice[:,cols[R2i]]
        R2iVec[R2iVec<0] = 0.0
        ax.plot(dataSlice[:,cols['dNdeta']]**(1.0/3.0), \
                R2iVec, \
                '-'+styles[KTIndex], \
                color=colors[KTIndex])
    ax.set_xlabel(r'$\left(dN_{\mathrm{ch}}/d\eta\right)^{1/3}$',
                    fontsize=16)
    ax.set_ylabel(R2iLabels[R2i], fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    #plt.show(block=False)
    #plt.show()
    fig.savefig("./"+R2i+"_vs_dNdeta.pdf")
    return None

def plot_R2i_vs_mult_dummy(data):
    fig, ax = plt.subplots()
    for KTIndex in range(data.shape[1]):
        dataSlice=data[:,KTIndex,:]
        ax.plot(dataSlice[:,cols['dNdeta']]**(1.0/3.0), \
                -1.0+0.0*dataSlice[:,cols['dNdeta']], \
                '-'+styles[KTIndex], \
                color=colors[KTIndex], \
                label=KTLabels[KTIndex])
    ax.set_ylim(bottom=0.0)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.legend(ncol=1, loc='center', fontsize=18)
    plt.tight_layout()
    fig.savefig("./R2i_vs_dNdeta_dummy.pdf")
    #plt.show(block=False)
    #plt.show()
    return None
    
    
def plot_CF_vs_KT(data, direction):
    fig, ax = plt.subplots()
    ax.axhline(1, color='k', linestyle='-', linewidth=1.0)
    for KTIndex in range(data.shape[2])[0:6]:
        dNdetaIndex = 0
        dataSlice=data[dNdetaIndex,CFcols[direction],KTIndex]
        ax.plot(dataSlice[:,2], \
                dataSlice[:,3], \
                '-'+styles[KTIndex], \
                color=colors[KTIndex])
        # include error bands
        plt.fill_between(dataSlice[:,2],\
                    dataSlice[:,3]-dataSlice[:,4],\
                    dataSlice[:,3]+dataSlice[:,4],\
                    alpha=0.25, edgecolor=colors[KTIndex],\
                    facecolor=colors[KTIndex])

    ax.set_xlabel(r'$q_{\mathrm{%(direction)s}}$ (GeV)' % {'direction': direction},
                    fontsize=16)
    ax.set_ylabel(r'$C$', fontsize=16)
    ax.set_xlim(left=-0.25, right=0.25)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    #plt.show()
    fig.savefig("./CF_vs_q" + direction + "_vs_KT.pdf")
    return None
    

def plot_CF_vs_dNdeta(data, direction):
    fig, ax = plt.subplots()
    ax.axhline(1, color='k', linestyle='-', linewidth=1.0)
    for dNdetaIndex in range(data.shape[0]):
        KTIndex = 0
        dataSlice=data[dNdetaIndex,CFcols[direction],KTIndex]
        ax.plot(dataSlice[:,2], \
                dataSlice[:,3], \
                '-'+styles[dNdetaIndex], \
                color=colors[dNdetaIndex])
        # include error bands
        plt.fill_between(dataSlice[:,2],\
                    dataSlice[:,3]-dataSlice[:,4],\
                    dataSlice[:,3]+dataSlice[:,4],\
                    alpha=0.25, edgecolor=colors[dNdetaIndex],\
                    facecolor=colors[dNdetaIndex])
    ax.set_xlabel(r'$q_{\mathrm{%(direction)s}}$ (GeV)' % {'direction': direction},
                    fontsize=16)
    ax.set_ylabel(r'$C$', fontsize=16)
    ax.set_xlim(left=-0.25, right=0.25)
    ax.set_ylim(top=1.4)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    #plt.show()
    fig.savefig("./CF_vs_q" + direction + "_vs_dNdeta.pdf")
    return None


if __name__ == "__main__":
    data = np.stack(tuple([get_dNchdeta_and_R2i_data(directory) for directory in sys.argv[1:]]))
    plot_R2i_vs_KT(data,'R2o')
    plot_R2i_vs_KT(data,'R2s')
    plot_R2i_vs_KT(data,'R2l')
    #plot_R2i_vs_KT_dummy(data)
    #plot_R2i_vs_mult(data,'R2o')
    #plot_R2i_vs_mult(data,'R2s')
    #plot_R2i_vs_mult(data,'R2l')
    #plot_R2i_vs_mult_dummy(data)
    data = np.stack(tuple([get_dNchdeta_and_CF_data(directory) for directory in sys.argv[1:]]))
    plot_CF_vs_KT(data, 'out')
    plot_CF_vs_KT(data, 'side')
    plot_CF_vs_KT(data, 'long')
    #plot_CF_vs_dNdeta(data, 'out')
    #plot_CF_vs_dNdeta(data, 'side')
    #plot_CF_vs_dNdeta(data, 'long')
    #pause()

# End of file
