#import matplotlib.pyplot as plt
import numpy as np
import sys

def get_event_class_mean_dNchdeta(filename):
    return np.mean(np.loadtxt(filename, usecols=4))

def get_directory_data(directory):
    mean_dNchdeta = get_event_class_mean_dNchdeta(directory + '/CF_results/event_class_multiplicities.dat')
    R2ij_vs_KT = np.loadtxt(directory + '/fit_results/HBTradii.dat')
    return np.c_[ np.broadcast_to( mean_dNchdeta, len(R2ij_vs_KT) ), R2ij_vs_KT ]


if __name__ == "__main__":
    data = np.stack(tuple([get_directory_data(directory) for directory in sys.argv[1:]]))
    print data.shape

# End of file