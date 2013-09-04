from __future__ import division 
import numpy as np
import mpmath as mp

if __name__ == "__main__":

    variance = []
    nvariance = []
    covariance = []
    correlation = []
    
    K = 250
    steps = 10000000
    i_in = 2
    
    for i in xrange(2,steps,10):
        #start i from 2
        v = np.exp(np.log(2/(2+3*K+K**2)) + (i-1)*np.log(K/(K+2))) - np.exp(-2*np.log(K) + 2*i*np.log(K/(K+1)))
        nv = np.sqrt(np.exp(np.log(2) + i*np.log(K/(K+2)) + (1-2*i)*np.log(K/(K+1))) - 1)
        #nv = np.exp(np.log(2/(K+1))+i*np.log((K+1)/(K+2))) - np.exp(i*np.log(K/(K+1))-np.log(K))
        variance.append([i,v])
        nvariance.append([i,nv])
    
    variance = np.array(variance)
    nvariance = np.array(nvariance)
    
    i = i_in
    vi = mp.exp(mp.log(2/(2+3*K+K**2)) + (i-1)*mp.log(K/(K+2))) - mp.exp(-2*mp.log(K) + 2*i*mp.log(K/(K+1)))
    for j in xrange(3,steps,10):
        #start j from 3    
        #this is wrong# c = mp.exp(-2*mp.log(K+1) + (i-1)*mp.log(K/(K+2)) + (j-1)*mp.log(K/(K+1))) - mp.exp(-2*mp.log(K) + (i+j)*mp.log(K/(K+1)))
        c = mp.exp((j-i)*mp.log(K)+(i-j-1)*mp.log(K+1)-i*mp.log(K+2)) - mp.exp(-2*mp.log(K) + (i+j)*mp.log(K/(K+1)))
        vj = mp.exp(mp.log(2/(2+3*K+K**2)) + (j-1)*mp.log(K/(K+2))) - mp.exp(-2*mp.log(K) + 2*j*mp.log(K/(K+1)))
        corr = c/mp.sqrt(vi*vj)
        covariance.append([j,c])
        correlation.append([j, corr])
    covariance = np.array(covariance)
    correlation = np.array(correlation)
    
    np.savetxt('variance_K{k}_steps{s}.dat'.format(k=K,s=steps),variance)
    np.savetxt('nvariance2_K{k}_steps{s}.dat'.format(k=K,s=steps),nvariance)
    np.savetxt('covariance_K{k}_steps{s}.dat'.format(k=K,s=steps),covariance)
    np.savetxt('correlation_K{k}_steps{s}.dat'.format(k=K,s=steps),correlation)