import matplotlib.pyplot as plt
import numpy as np
import sys

columnLabels = ['dNdeta','KT',\
                'lambda','R2o','R2s','R2l',\
                'err(lambda)','err(R2o)','err(R2s)','err(R2l)']
multLabels = [r'$N_{ch}=1-11$', r'$N_{ch}=12-16$', r'$N_{ch}=17-22$', \
              r'$N_{ch}=23-28$', r'$N_{ch}=29-34$', r'$N_{ch}=35-41$', \
              r'$N_{ch}=42-51$', r'$N_{ch}=52-151$']

cols = dict(zip(columnLabels,range(len(columnLabels))))
colors = ['blue','red','green','purple','cyan','magenta','orange','yellow']
styles = ['o','*','v','^','s','D','P','X']


def get_event_class_mean_dNchdeta(filename):
    return np.mean(np.loadtxt(filename, usecols=4))

def get_directory_data(directory):
    mean_dNchdeta = get_event_class_mean_dNchdeta(directory + '/CF_results/event_class_multiplicities.dat')
    R2i_vs_KT = np.loadtxt(directory + '/fit_results/HBTradii.dat', \
                           usecols=(0,3,4,5,6,10,11,12,13))
    return np.c_[ np.broadcast_to( mean_dNchdeta, len(R2i_vs_KT) ), R2i_vs_KT ]

def plot_R2i_vs_KT(data, R2i):
    fig, ax = plt.subplots()
    for multIndex, dataSlice in np.ndenumerate(data):
        ax.plot(dataSlice[:,cols['KT']], \
                dataSlice[:,cols[R2i]], \
                '-'+styles[multIndex], color=colors[multIndex])
    fig.savefig("./"+R2i+"_vs_KT.pdf")
    return None

def plot_R2i_vs_mult(data, R2i):
    dataT = np.swapaxes(data,0,1)
    fig, ax = plt.subplots()
    for KTIndex, dataTSlice in np.ndenumerate(dataT):
        ax.plot(dataTSlice[:,cols['dNdeta']]**(1.0/3.0), \
                dataTSlice[:,cols[R2i]], \
                '-'+styles[KTIndex], color=colors[KTIndex])
    fig.savefig("./"+R2i+"_vs_dNdeta.pdf")
    return None

if __name__ == "__main__":
    data = np.stack(tuple([get_directory_data(directory) for directory in sys.argv[1:]]))
    plot_R2i_vs_KT(data,'R2o')
    plot_R2i_vs_KT(data,'R2s')
    plot_R2i_vs_KT(data,'R2l')
    plot_R2i_vs_mult(data,'R2o')
    plot_R2i_vs_mult(data,'R2s')
    plot_R2i_vs_mult(data,'R2l')
    #print(data.shape)

# End of file
