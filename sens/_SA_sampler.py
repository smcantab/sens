import numpy as np
from scipy.special import gammaln

from pele.utils.rotations import vec_random_ndim

from sens.src.weighted_pick import weighted_pick_cython


#def vector_random_uniform_hypersphere(k):
#    """return a vector sampled uniformly in a hypersphere of dimension k"""
#    u = vec_random_ndim(k)
#    #draw the magnitude of the vector from a power law density:
#    #draws samples in [0, 1] from a power distribution with positive exponent k - 1.
#    p = np.random.power(k)
#    return p * u
#
#def sample_uniformly_in_basin_harmonic(m, Emax, k):
#    """ Returns a configuration with energy less than Emax sampled uniformly from m-th basin    
#        this assumes the harmonic approximation and is exact in the harmonic approximation
#        
#        Parameters
#        ----------
#        
#        m : Minimum object
#            normal mode frequencies and associated eigenvectors
#        Emax : float
#            energy upper bound
#        k : integer
#            number of degrees of freedom
#    """
#    nm = m.hessian_eigs[0]
#    evals = nm.eigenvalues
#    vectors = nm.eigenvectors
#    nzero = len(evals) - k
#
#    # get uniform random k dimensional unit vector
#    f = vec_random_ndim(k)
#
#    # the target increase in energy is sampled from a power law distribution
#    dE_target = np.random.power(k-1) * (Emax - m.energy) ####CHANGED TO K-1 
#    
#    # scale f according to dE_target
#    f *= np.sqrt(2. * dE_target)
#
#    # create the random displacement vector
#    dx = np.zeros(m.coords.shape)
#    for i in range(k):
#        if evals[i+nzero] > 1e-4:
#            dx += f[i] * vectors[:,i+nzero] / np.sqrt(evals[i+nzero])
#    
#    return m.coords + dx
#    
#
#def sample_uniformly_in_basin(m, Emax, potential, k):
#    """return a configuration with energy less than Emax sampled uniformly from the basin defined by m
#    
#        Parameters
#        ----------
#        
#        m : Minimum object
#            normal mode frequencies and associated eigenvectors
#        Emax : float
#            energy upper bound
#        Potential: type Potential
#            potential energy function, takes coordinates and returns energy configuration
#        k : integer
#            number of degrees of freedom
#    """
#    # displace randomly from the minimum according to the eigenvalues and eigenvectors
#    E = Emax + 1.
#    count = 0
#    while E > Emax:
#        coords = sample_uniformly_in_basin_harmonic(m, Emax, k)
#        E = potential.getEnergy(coords)
#        
#        # print some stuff
#        stepsize = np.linalg.norm(coords - m.coords)
#        print "created structure with energy", E, "Emax", Emax, "Emin", m.energy, count, stepsize
#        count += 1
#        
#    
#    # now do a short monte carlo sampling to improve unbiased sampling
#    
#    return coords, E
#
#def compute_log_phase_space_volume(m, Emax, k):
#    """return the log (base e) of the phase space volume of minima m up to energy Emax
#    
#    Notes
#    -----
#    
#    V = Integral from m.energy to Emax of the harmonic density of states DoS(E)
#    
#        DoS(E) = 2*N!*(E - m.energy)^(k-1) / (Gamma(k) * prod_freq * O_k)
#    
#    After integration between m.energy to Emax one finds the volume:  
#        
#    V = 2*N!*(Emax - m.energy)^k / (np.gamma(k+1) * prod_freq * O_k)
#    
#    note: gammaln(N+1) approximates the ln(N!)
#    
#    Parameters
#    ----------
#    
#    k : integer
#        number of vibrational degrees of freedom
#    Gamma : scipy.special function 
#        gamma function
#    prod_freq : float
#        the product of the frequencies (from the eigenvalues of the Hessian)
#    O_k : integer
#        the order of the symmetry point group
#    """
#    from numpy import log
#    
#    if m.energy > Emax:
#        raise ValueError("Emax (%g) must be greater than the energy of the minimum (%g)" % (Emax, m.energy))
#    logV = log(2) + gammaln(N+1) + k * log(Emax - m.energy) - gammaln(k+1) - m.fvib/2. - log(m.pgorder)
#    return logV
#
#
#def sample_minimum(minima, Emax, k):
#    """return a minimum sampled uniformly with weight according to phase space volume
#    
#    Parameters
#    ----------
#    minima : list of Mimumum objects
#    Emax : float
#        the maximum energy for the phase space volume calculation
#    k : int
#        the number of degrees of vibrational freedom (3*N-6 for atomic clusters)
#    """
#    # calculate the harmonic phase space volume of each minima and store it in list `weights`
#    lweights = []
#    minima2 = []
#    for m in minima:
#        if m.energy < Emax:
#            lV = compute_log_phase_space_volume(m, Emax, k)
#            lweights.append(lV)
#            minima2.append(m)
#    lweights = np.array(lweights)
#    weights = np.exp(lweights - np.max(lweights))
#    
#    # select a minimum uniformly given `weights`
#    # print "weights", weights[:10]
#    index = weighted_pick_cython(weights)
#    # print index, len(weights), len(minima)
#    m = minima2[index]
#    return m
#
#def sample_from_database(system, minima, Emax):
#    """return a configuration sampled uniformly from a database of minima according to the harmonic approximation up to energy Emax
#    
#    Parameters
#    ----------
#    
#    system : pele System
#        is the particular system of interest, say LJCluster
#    minima: array of Minimum objects
#    Emax : float
#        energy upper bound
#    
#    """
#    m = sample_minimum(minima, Emax, system.k)
#    
#    # sample configuration uniformly from the basin of minima m 
#    coords, E = sample_uniformly_in_basin(m, Emax, system.get_potential(), system.k)
#
#    return coords, E
 
#def populate_database(system, db, niter=1000):
#    """return a database with all important low energy minima 
#    
#    """
#    # use basinhopping to find the low energy minima and store them in a database
#    bh = system.get_basinhopping(database=db)
#    bh.run(niter)


class SASampler(object):
    """this class will manage the sampling of configurations from a database of minima
    
    in particular it will precompute values so they need not be recalculated every time
    """
    def __init__(self, minima, k):
        self.minima = minima
        self.k = k
        self.gammalnk = gammaln(self.k)
        self.lVol_prefactor = self.precompute_log_phase_space_volume_prefactor()
        
        self.energyvec = np.array([m.energy for m in self.minima])
        self.lVol_prefactor_vec = np.array([self.lVol_prefactor[m] for m in self.minima])
        
    def log_phase_space_volume_prefactor(self, m):
        """return the log of the part of the volume excluding Emax
        
        Notes
        -----
        
        fvib: log product of squared frequencies
        
        pgorder: point group order
            
        """
        #return - np.log(self.k) - self.gammalnk - m.fvib/2. - np.log(m.pgorder)
        return - m.fvib/2. - np.log(m.pgorder)
    
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
        nm = m.hessian_eigs[0]
        evals = nm.eigenvalues
        vectors = nm.eigenvectors
        k = self.k
        nzero = len(evals) - k
    
        # get uniform random k dimensional unit vector
        f = vec_random_ndim(k)
    
        # the target increase in energy is sampled from a power law distribution
        dE_target = np.random.power(k-1) * (Emax - m.energy)
        
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