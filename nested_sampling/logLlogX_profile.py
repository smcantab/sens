import argparse
import numpy as np
import copy
from itertools import cycle
import matplotlib
from matplotlib import rc
import os.path
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import mpmath as mp

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute Cv stdev")
    parser.add_argument("K", type=int, help="number of replicas")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with energies")
    parser.add_argument("-P", type=int, help="number of cores for parallel run", default=1)
    parser.add_argument("--live_not_stored", action="store_true", help="turn this flag on if you're using a set of data that does not contain the live replica.",default=False)
    args = parser.parse_args()
    
    fname = args.fname
    K = args.K
    P = args.P
        
    all_data = [[] for i in xrange(len(fname))]
    all_labels = np.array([os.path.splitext(name)[0] for name in fname])
    for name,i in zip(fname,xrange(len(fname))):
        data = np.genfromtxt(name)
        all_data[i] = data
        
    all_data = np.array(all_data)
    
    if args.live_not_stored is True:
        c = 0
    else:
        c = K   
    
    for data, name in zip(all_data, all_labels):
        tmp = np.zeros((len(data),2))
        tmp[:,1] = data
        lX = 0
        L = np.shape(tmp)[0]
        
        for i in xrange(L-c):
            lX += np.log(K-i%P) - np.log(K+1-i%P)
            tmp[i,0] = np.exp(lX)
        
        if args.live_not_stored is False:
            for i in xrange(K):
                lX += np.log(K-i) - np.log(K+1-i)
                tmp[L-c+i,0] = np.exp(lX)
        
        np.savetxt('{n}_logLlogX.dat'.format(n=name),tmp)
    
    
                
            
        