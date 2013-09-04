from __future__ import division
import numpy as np
import matplotlib
from matplotlib import rc
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy import special
from itertools import cycle

#######################SET LATEX OPTIONS###################
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})
#rc('text.latex',preamble=r'\usepackage{times}')
plt.rcParams.update({'font.size': 14})
plt.rcParams['xtick.major.pad'] = 8
plt.rcParams['ytick.major.pad'] = 8
##########################################################

def Ei(i,k,K,P,Eg,V):
    alpha = 1-P/(K+1)
    if i is 0:
        return Eg
    else:
        return (Eg-V)*alpha**(2*i/k) + V*(alpha**(2*i/k)-1)/(alpha**(2/k)-1)

def dPb(P,i,k,K,Eg,Va,Vb,va,vb,oa,ob):
    E = Ei(i,k,K,P,Eg,Va)
    Pb = (E-Vb)**(k/2) * (va**k)*oa/((E-Vb)**(k/2) * (va**k)*oa + (E-Va)**(k/2) * (vb**k)*ob)
    return Pb

def n_best(k,K,P,V,Eg,v,d):
    alpha = 1-P/(K+1)
    G = ((Eg-V)**(k/2)) / ((2*np.pi)**(k/2) * special.gamma(k/2 +1) * v)
    n = (np.log(G) + d*np.log(10))/np.log(1/alpha)
    return int(np.rint(n))

def Pb(P,i,k,K,Eg,Va,Vb,va,vb,oa,ob):
    nb = n_best(k,K,P,Va,Eg,va,d)
    sumPb = 0
    for i in xrange(nb):
        sumPb += dPb(P,i,k,K,Eg,Va,Vb,va,vb,oa,ob)
    return P*sumPb

if __name__ == "__main__":
    
    d = 9 #n_best to get within 10^-9  
    k = 3
    #K = 1000
    P = 1   
    Eg = 500
    Va = 100
    Vb = 0
    va = 1
    vb = 10
    oa = 1
    ob = oa
    
    Ga = ((Eg-Va)**(k/2)) / ((2*np.pi)**(k/2) * special.gamma(k/2 +1) * va)
    Gb = ((Eg-Vb)**(k/2)) / ((2*np.pi)**(k/2) * special.gamma(k/2 +1) * vb)
    assert Ga > Gb
    print "Gsub",Ga,"Gdom",Gb 
    
    #nb = n_best(k,K,P,Va,Eg,va,d)
    #print "nb", nb
    #sumPb = Pb(nb,P,i,k,K,Eg,Va,Vb,va,vb,oa,ob) 
    #print "Pds",1/sumPb,"a",K/sumPb,"nb",nb
    #####################LINE STYLE CYCLER####################
    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)
    cm = plt.get_cmap('Dark2')
    ##########################################################
    
    all_data = 4
    Kmin = 10
    Kmax = 1000
    nK = (Kmax - Kmin) 
    dK = int((Kmax - Kmin)/nK)
    karray = np.zeros(nK)
    parray = np.zeros(nK)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_color_cycle([cm(1.*i/all_data) for i in xrange(all_data)])
    P = 2 # remove this if want to run this for Va
    for i in xrange(all_data):
        K = Kmin
        for i in xrange(nK):
            karray[i] = K
            parray[i] = -K*np.log(1-0.999)/Pb(P,i,k,K,Eg,Va,Vb,va,vb,oa,ob)
            K += dK
        ax.plot(karray, parray, next(linecycler),label='$\mathcal{P}=%s$'%(P), linewidth=1.8)
        P += 2
        #ax.plot(karray, parray, next(linecycler), label=r'$V_a={V}$'.format(V=Va), linewidth=1.8)
        #Va -= 25
        print P
    #ax.set_title(r'$\kappa = 3$ $E_g = 500$ $V_a = 100$ $V_b = 0$ $\overline{v}_a=1$ $\overline{v}_b=10$ $o_b = o_a = 1$ $d=9$ $\phi=0.999$')
    ax.set_title('(b)')
    ax.set_xlabel('$K$')
    ax.set_ylabel('$a = KP_{DS}$')
    ax.legend(frameon=False,loc=4)
    plt.show()
    fig.savefig('ideal_pds_P.eps')
    
