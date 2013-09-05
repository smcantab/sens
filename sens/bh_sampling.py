"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np
import random
from scipy.special import gamma, gammaln
from scipy.misc import factorial


from pele.utils.rotations import vec_random_ndim
from pele.utils.hessian import sort_eigs, get_eig
from pele.thermodynamics import logproduct_freq2, normalmodes

from nested_sampling import NestedSampling, Replica

from sens.database_eigenvecs import HessianEigs
from sens._SA_sampler import SASampler




class ConfigTestError(StandardError):
    pass

class NestedSamplingSA(NestedSampling):
    """overload get_starting_configuration() in order to introduce sampling from known minima
    
    Parameters
    ----------
    system : pele system object
    nreplicas : int
        number of replicas
    takestep : callable
        object to do the step taking.  must be callable and have attribute takestep.stepsize
    minima : list of Minimum objects
    """
    def __init__(self, system, nreplicas, mc_runner, minima, 
                  minprob=None, energy_offset=None,
                  **kwargs):
        super(NestedSamplingSA, self).__init__(system, nreplicas, mc_runner, **kwargs)
        self.minima = minima
        self.bh_sampler = SASampler(self.minima, self.system.k)
        if minprob is None:
            raise ValueError("minprob cannot be None")
        self.minprob = minprob
        if energy_offset is None:
            self.energy_offset = 2.5
        else:
            self.energy_offset = energy_offset
        
        self.count_sampled_minima = 0
        
    def get_starting_configuration_minima_HA(self, Emax):
        """using the Harmonic Approximation sample a new configuration starting from a minimum sampled uniformly according to phase space volume
        
        Notes
        -----
        displace minimum configuration along its eigenvectors
        
        """
        m = self.bh_sampler.sample_minimum(Emax)        
        x = self.bh_sampler.sample_coords_from_basin(m, Emax)
        pot = self.system.get_potential()
        e = pot.getEnergy(x)
        return x, e
    
    def get_starting_configuration_minima_single(self, Emax):
        """return a minimum sampled uniformly according to phase space volume"""
        m = self.bh_sampler.sample_minimum(Emax)
        x, e = m.coords, m.energy
        self.system.center_coords(x)
        if True:
            accept_tests = self.system.get_config_tests()
            for test in accept_tests:
                t = test(coords=x)
                if not t:
                    print "warning: minimum from database failed configuration test", m._id, m.energy
                    raise ConfigTestError()
        return x, e

    def get_starting_configuration_minima(self, Emax):
        """return a minimum sampled uniformly according to phase space volume

        Notes
        -----
        some of the minima in the database don't satisfy the configuration
        checks.  So we sample over and over again until we get one that
        passes

        """
        self.count_sampled_minima += 1
        while True:
            try:
                return self.get_starting_configuration_minima_single(Emax)
            except ConfigTestError:
                pass

    def onset_prob_func(self, Emax):
        """return the probability of sampling from a minimum
        
        The probability depends on Emax. For high Emax the probability is 0.  the probability
        increases as Emax get's lower, reaching a maximum of self.minprob.  

        value of Emax where the probabilty starts to get large is energy_max_database + energy_offset where
        energy_max_database is the maximum energy minimum in the database.  This behavior
        can be adjusted with parameter energy_offset.
        
        parameter b determines the speed at which it turns on.
        
        """
        if not hasattr(self, "_energy_max_database"):
            self._energy_max_database = float(max([m.energy for m in self.minima]))
        max_prob = float(self.minprob)
        energy_onset_width = 1.
        dE = self._energy_max_database + self.energy_offset - Emax
        f = dE / energy_onset_width
        if f > 100:
            onset_prob = 0.
        else:
            onset_prob = max_prob / ( 1. + np.exp(-f))
        return onset_prob
    
    def get_starting_configurations(self, Emax):
        """this function overloads the function in NestedSampling"""
        # choose a replica randomly
        configs = self.get_starting_configurations_from_replicas()
        # replace each starting configuration with a one chosen
        # from the minima with probability prob
        onset_prob = self.onset_prob_func(Emax)
        prob = onset_prob / float(self.nreplicas)
        for i in range(len(configs)):
            if np.random.uniform(0,1) < prob:
                x, energy = self.get_starting_configuration_minima(Emax)
                configs[i] = Replica(x, energy, from_random=False)
                if self.verbose:
                    print "sampling from minima, E minimum:", energy, "with probability:", prob
        return configs

if __name__ == "__main__":
    # define the system
    from lj_run import LJClusterNew
    natoms = 13
    system = LJClusterNew(natoms)

    db = system.create_database("lj%d.db" % (natoms))
#    if True:
#        populate_database(system, db, niter=100)
    
    print "pgorder", db.minima()[0].pgorder
    print "fvib", db.minima()[0].fvib
    get_thermodynamic_information(system, db)
    print "pgorder", db.minima()[0].pgorder
    print "fvib", db.minima()[0].fvib
    
    Emin = db.minima()[0].energy
    Emax = Emin + 1.

    k = system.k
    for i in range(10):
        coords, E = sample_from_database(system, db, Emax)
    
