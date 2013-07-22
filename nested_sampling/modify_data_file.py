import argparse
import numpy as np
import copy
from itertools import cycle
import matplotlib
from matplotlib import rc
import os.path
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

##should write it in such a way that the operations are carried out with priority determined by the order in which
##they are listed in the command line, rather than set priority

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute Cv stdev", 
                                     epilog="must write file name followed by a label string, otherwise raise the flag --nolabels")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with heat capacity followed by label")
    parser.add_argument("--column", type=int, help="set data column to modify",default=None)
    parser.add_argument("-m", type=float, help="factor by which multiply whole column, priority=1",default=1.)
    parser.add_argument("-d", type=float, help="factor by which divide whole column, priority=1",default=1.)
    parser.add_argument("-a", type=float, help="factor by which add whole column, priority=1",default=0.)
    parser.add_argument("-s", type=float, help="factor by which subtract whole column, priority=1",default=0.)
    parser.add_argument("--di", type=float, help="divide this quantity by the values of the selected column, priority=2",default=0)
    parser.add_argument("--log", action="store_true", help="take natural log of selected column, priority=3",default=False)
    parser.add_argument("--exp", action="store_true", help="take exp of selected column, priority=4",default=False)
    args = parser.parse_args()
    
    fname = args.fname
    print fname
    
    column = args.column
    add = args.a
    multiply = args.m
    divide = args.d
    subtract = args.s
    inv_divide = args.di
    flog = args.log
    fexp = args.exp
    
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
    if inv_divide is not None:
        data[:][column] = inv_divide / data[:][column]
    if flog is True:
        data[:][column] = np.log(data[:][column])
    if fexp is True:
        data[:][column] = np.exp(data[:][column])
 
    all_data = np.array(all_data)
    
    if args.column is not None:
        for data in all_data:
            data[:][column] = (data[:][column] + add - subtract) * multiply / divide
    
    for data, name in zip(all_data, all_labels):
        np.savetxt('{n}_modified.dat'.format(n=name),data)
   
