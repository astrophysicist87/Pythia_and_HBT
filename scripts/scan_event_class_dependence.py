import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

columnLabels = ['dNdeta','KT',\
                'lambda','R2o','R2s','R2l',\
                'err(lambda)','err(R2o)','err(R2s)','err(R2l)']
multLabels = [r'$N_{ch}=1-11$', r'$N_{ch}=12-16$', r'$N_{ch}=17-22$', \
              r'$N_{ch}=23-28$', r'$N_{ch}=29-34$', r'$N_{ch}=35-41$', \
              r'$N_{ch}=42-51$', r'$N_{ch}=52-151$']
KTLabels = [r'$K_T=0-100$ MeV', r'$K_T=100-200$ MeV', r'$K_T=200-300$ MeV', \
            r'$K_T=300-400$ MeV', r'$K_T=400-500$ MeV', r'$K_T=500-600$ MeV', \
            r'$K_T=600-700$ MeV', r'$K_T=700-800$ MeV']

R2iLabels = dict({'R2o':r'$R^2_o$ fm$^2$','R2s':r'$R^2_s$ fm$^2$','R2l':r'$R^2_l$ fm$^2$'})

cols = dict(zip(columnLabels,range(len(columnLabels))))
colors = ['blue','red','green','purple','cyan','magenta','orange','yellow']
styles = ['o','*','v','^','s','D','P','X']

def pause():
    programPause = raw_input("")

def get_event_class_mean_dNchdeta(filename):
    return np.mean(np.loadtxt(filename, usecols=4))

def get_directory_data(directory):
    mean_dNchdeta = get_event_class_mean_dNchdeta(directory + '/CF_results/event_class_multiplicities.dat')
    R2i_vs_KT = np.loadtxt(directory + '/fit_results/HBTradii.dat', \
                           usecols=(0,3,4,5,6,10,11,12,13))
    return np.c_[ np.broadcast_to( mean_dNchdeta, len(R2i_vs_KT) ), R2i_vs_KT ]

def plot_R2i_vs_KT(data, R2i):
    fig, ax = plt.subplots()
    for multIndex in range(len(data)):
        dataSlice=data[multIndex]
        ax.plot(dataSlice[:,cols['KT']], \
                dataSlice[:,cols[R2i]], \
                '-'+styles[multIndex], \
                color=colors[multIndex])
    #plt.show(block=False)
    #plt.show()
    ax.set_xlabel(r'$K_T$ (GeV)')
    ax.set_ylabel(R2iLabels[R2i])
    fig.savefig("./"+R2i+"_vs_KT.pdf")
    return None

def plot_R2i_vs_KT_dummy(data):
    fig, ax = plt.subplots()
    for multIndex in range(len(data)):
        dataSlice=data[multIndex]
        ax.plot(dataSlice[:,cols['KT']], \
                -1.0+0.0*dataSlice[:,cols['KT']], \
                '-'+styles[multIndex], \
                color=colors[multIndex], \
                label=multLabels[multIndex])
    ax.set_ylim(bottom=0.0)
    ax.legend(ncol=2,loc='best')
    fig.savefig("./R2i_vs_KT_dummy.pdf")
    #plt.show(block=False)
    #plt.show()
    return None

def plot_R2i_vs_mult(data, R2i):
    fig, ax = plt.subplots()
    for KTIndex in range(data.shape[1]):
        dataSlice=data[:,KTIndex,:]
        ax.plot(dataSlice[:,cols['dNdeta']]**(1.0/3.0), \
                dataSlice[:,cols[R2i]], \
                '-'+styles[KTIndex], \
                color=colors[KTIndex])
    #plt.show(block=False)
    #plt.show()
    ax.set_xlabel(r'$\left(dN_{\mathrm{ch}}/d\eta\right)^{1/3}$')
    ax.set_ylabel(R2iLabels[R2i])
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
    ax.legend(ncol=2,loc='best')
    fig.savefig("./R2i_vs_dNdeta_dummy.pdf")
    #plt.show(block=False)
    #plt.show()
    return None

if __name__ == "__main__":
    data = np.stack(tuple([get_directory_data(directory) for directory in sys.argv[1:]]))
    plot_R2i_vs_KT(data,'R2o')
    plot_R2i_vs_KT(data,'R2s')
    plot_R2i_vs_KT(data,'R2l')
    plot_R2i_vs_KT_dummy(data)
    plot_R2i_vs_mult(data,'R2o')
    plot_R2i_vs_mult(data,'R2s')
    plot_R2i_vs_mult(data,'R2l')
    plot_R2i_vs_mult_dummy(data)
    #print(data.shape)
    #pause()

# End of file
