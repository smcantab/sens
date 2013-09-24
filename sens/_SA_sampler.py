import numpy as np
from scipy.special import gammaln

from pele.utils.rotations import vec_random_ndim
from nested_sampling.utils.rotations import vector_random_uniform_hypersphere

from sens.src.weighted_pick import weighted_pick_cython


class _UnboundNormalModes(object):
    def __init__(self, nm):
        self.freqs = nm.freqs.copy()
        self.vectors = nm.vectors.copy()

class _UnboundMinimum(object):
    """represent a minimum object unbound from the database"""
    def __init__(self, m):
#        self._id = m._id
        self.energy = m.energy
        self.coords = m.coords.copy()
        self.fvib = m.fvib
        self.pgorder = m.pgorder
        try:
            if m.normal_modes is not None:
                self.normal_modes = _UnboundNormalModes(m.normal_modes)
            else:
                self.normal_modes = None
        except AttributeError:
            self.normal_modes = None
    
#    def __eq__(self, m2):
#        assert self._id is not None
#        assert m2._id is not None
#        return self._id == m2._id
#
#    def __hash__(self):
#        assert self._id is not None
#        return self._id


class SASampler(object):
    """this class will manage the sampling of configurations from a database of minima
    
    It will precompute values so they need not be recalculated every time
    """
    def __init__(self, minima, k, copy_minima=True, center_minima=False):
        if copy_minima or center_minima:
            self.minima = [_UnboundMinimum(m) for m in minima]
            if center_minima:
                # js850: this is a quick hack.  this should be done more elegantly
                for m in self.minima:
                    x = m.coords.reshape([-1,3])
                    com = x.sum(0) / x.shape[0]
                    x = x - com[np.newaxis,:]
                    x = x.reshape(-1)
                    m.coords[:] = x[:] 
        else:
            self.minima = minima
        self.minima = sorted(self.minima, key=lambda m:m.energy)

        self.k = k
        self.gammalnk = gammaln(self.k)
        self.lVol_prefactor = self.precompute_log_phase_space_volume_prefactor()
        
        self.energyvec = np.array([m.energy for m in self.minima])
        self.lVol_prefactor_vec = np.array([self.lVol_prefactor[m] for m in self.minima])
        
    def log_phase_space_volume_prefactor(self, m):
        """return the log of the part of the volume excluding Emax
        
        Notes
        -----
        
        fvib: log product of *squared* frequencies
        
        pgorder: point group order
            
        """
        #return - np.log(self.k) - self.gammalnk - m.fvib/2. - np.log(m.pgorder)
        return - 0.5 * m.fvib - np.log(m.pgorder)
    
    def precompute_log_phase_space_volume_prefactor(self):
        """return a dictionary of the log phase space volume prefactors
        
        precompute the parts of the log phase space volume that don't depend on Emax
        """ 
        return dict([(m, self.log_phase_space_volume_prefactor(m)) for m in self.minima]) 

    def log_phase_space_volume(self, m, Emax):
        """return the log phase space volume of minimum m up to energy Emax
        
        Notes
        -----
        Only the terms that depend on which minimum is chosen are included.  The
        other terms would be canceled out upon normalization
        
        V = Integral from m.energy to Emax of the harmonic density of states DoS(E)
        
            DoS(E) = 2*N!*(E - m.energy)^(k-1) / (Gamma(k) * prod_freq * O_k)
        
        After integration between m.energy to Emax one finds the volume:  
            
        V = 2*N!*(Emax - m.energy)^k / (np.gamma(k+1) * prod_freq * O_k)
        
        note: all constant multiplicative factors are not necessary as they cancel out when renormalising the weights,
            thus reducing the V to:
        V = (Emax - m.energy)^k / (prod_freq * O_k)
               
        """
        return (self.k / 2)  * np.log(Emax - m.energy) + self.lVol_prefactor[m]

    def sample_coords_from_basin(self, m, Emax):
        """Returns a configuration with energy less than Emax sampled uniformly from the basin of a minimum
        
        this assumes the harmonic approximation and is exact in the harmonic approximation
        
        Parameters
        ----------
        m : Minimum object
            normal mode frequencies and associated eigenvectors
        Emax : float
            energy upper bound
        k : integer
            number of degrees of freedom
        
        Notes
        -----
        in real system, even ones that fit quite well to the harmonic approximation this very often generates 
        configurations with energy greater than Emax.  
        """
        nm = m.normal_modes
        evals = nm.freqs
        vectors = nm.vectors
        k = self.k
        nzero = len(evals) - k
    
        # get uniform random k dimensional unit vector
        f = vector_random_uniform_hypersphere(k)
    
        # the target increase in energy is sampled from a power law distribution
        # js850 sep 13> dE_target = np.random.power(k-1) * (Emax - m.energy)
        dE_target = (Emax - m.energy)
        
        # scale f according to dE_target
        f *= np.sqrt(2. * dE_target) #TODO check prefactor
    
        # create the random displacement vector
        dx = np.zeros(m.coords.shape)
        for i in range(k):
            if evals[i+nzero] > 1e-4:
                dx += f[i] * vectors[:,i+nzero] / np.sqrt(evals[i+nzero])
        
        return m.coords + dx
    
    def _compute_weights(self, Emax):
        """compute weights of minima from phase space volumes up to energy Emax  
        
        Notes
        -----
        Only the terms that depend on which minimum is chosen are included.  The
        other terms would be canceled out upon normalization

        This version does all calculations with numpy vectors to improve speed
        
        V = Integral from m.energy to Emax of the harmonic density of states DoS(E)
        
            DoS(E) = 2*N!*(E - m.energy)^(k-1) / (Gamma(k) * prod_freq * O_k)
        
        After integration between m.energy to Emax one finds the volume:  
            
        V = 2*N!*(Emax - m.energy)^k / (np.gamma(k+1) * prod_freq * O_k)
        
        note: all constant multiplicative factors are not necessary as they cancel out when renormalising the weights,
            thus reducing the V to:
        V = (Emax - m.energy)^k / (prod_freq * O_k)
        
        """
        i_max = np.searchsorted(self.energyvec, Emax, side='right')
        #indices = np.where(self.energyvec < Emax)[0]
        #minima2 = [self.minima[i] for i in indices]
        minima2 = self.minima[:i_max]
        lweights = (float(self.k)/2.) * np.log(Emax - self.energyvec[:i_max]) + self.lVol_prefactor_vec[:i_max]
#        print "energyvec", self.energyvec
#        print "lvol_pref", self.lVol_prefactor_vec
        weights = np.exp(lweights - np.max(lweights))
        return minima2, weights
    
    def compute_weights(self, Emax):
        if not hasattr(self, "_weights_Emax"):
            self._weights_Emax = None
        if self._weights_Emax == Emax:
            return self._weights_minima, self._weights
        else:
            self._weights_Emax = Emax
            self._weights_minima, self._weights = self._compute_weights(Emax)
            return self._weights_minima, self._weights
            
            
    def compute_weights_slow(self, Emax):
        """compute weights from phase space volumes
        
        This version is slow because it uses attribute and dictionary lookups
        """
        minima2 = [m for m in self.minima if m.energy < Emax]
        lweights = [self.log_phase_space_volume(m, Emax) for m in minima2]
        lweights = np.array(lweights)
        weights = np.exp(lweights - np.max(lweights))
        return minima2, weights


    def sample_minimum(self, Emax):
        """return a minimum sampled uniformly with weight according to phase space volume
        
        Parameters
        ----------
        Emax : float
            the maximum energy for the phase space volume calculation
        """
        # calculate the harmonic phase space volume of each minima and store it in list `weights`
        minima2, weights = self.compute_weights(Emax)
        
        # select a minimum uniformly given `weights`
        index = weighted_pick_cython(weights)
        m = minima2[index]
        return m
    
    def sample_coords(self, Emax):
        m = self.sample_minimum(Emax)
        coords = self.sample_coords_from_basin(m, Emax)
        return m, coords
    
    def compute_energy(self, x, m, x0=None):
        """compute the harmonic energy of configuration x with respect to minimum m
        
        Notes
        -----
        x must be aligned with m.coords.  This means that trivial degrees of freedom like
        translation and rotation must be accounted for because they are not symmetries in the HSA.
        Depending on which degrees of freedom there are this alignment can be non-trivial, but the tools in
        pele.mindist can help.   
        """
        if x0 is None:
            x0 = m.coords
        dx = x - x0
        nm = m.normal_modes
        freqs = nm.freqs
        vectors = nm.vectors
        nzero = freqs.size - self.k
        
        energy = sum([0.5 * freqs[i] * np.dot(dx, vectors[:,i])**2 
                      for i in xrange(nzero, freqs.size)])
        
        energy += m.energy
        
        return energy
            
            
