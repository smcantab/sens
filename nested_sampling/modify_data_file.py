import argparse
import numpy as np
import copy
from itertools import cycle
import matplotlib
from matplotlib import rc
import os.path
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute Cv stdev", 
                                     epilog="must write file name followed by a label string, otherwise raise the flag --nolabels")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with heat capacity followed by label")
    parser.add_argument("--column", type=int, help="set data column to modify",default=None)
    parser.add_argument("-m", type=float, help="factor by which multiply whole column",default=1.)
    parser.add_argument("-d", type=float, help="factor by which divide whole column",default=1.)
    parser.add_argument("-a", type=float, help="factor by which add whole column",default=0.)
    parser.add_argument("-s", type=float, help="factor by which subtract whole column",default=0.)
    args = parser.parse_args()
    
    fname = args.fname
    
    column = args.column
    add = args.a
    multiply = args.m
    divide = args.d
    subtract = args.s
    
    
    ####################################################################################################
    #DEAL WITH INPUT
        
    all_data = [[] for i in xrange(len(fname))]
    all_labels = np.array([os.path.splitext(name)[0] for name in fname])
    for name,i in zip(fname,xrange(len(fname))):
        data = np.genfromtxt(name)
        #dshape = np.shape(data)
        #data = np.hstack((data, np.zeros((data.shape[0], 1), dtype=data.dtype)))
        data = np.array(data)
        all_data[i] = data
            
    all_data = np.array(all_data)
    
    if args.column is not None:
        for data in all_data:
            data[:][column] = (data[:][column] + add - subtract) * multiply / divide
    
    for data, name in zip(all_data, all_labels):
        np.savetxt('{n}_modified.dat'.format(n=name),data)
    