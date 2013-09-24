"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np

from nested_sampling import NestedSampling, Replica

from sens import SASampler




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
                  minprob=None, energy_offset=None, copy_minima=True, 
                  center_minima=False, **kwargs):
        super(NestedSamplingSA, self).__init__(system, nreplicas, mc_runner, **kwargs)
        self.minima = minima
        self.nreplicas = nreplicas
        self.bh_sampler = SASampler(self.minima, self.system.k, copy_minima, center_minima)
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
        energy_max_database is the maximum energy minimum in the database.  This behaviour
        can be adjusted with parameter energy_offset.
        
        parameter <energy_onset_width> determines the speed at which it turns on.
        
        the optimal probability of sampling from the database scales approximately as 1/nreplicas, 
        this has been shown analytically (see Stefano's thesis pp.29-32). 
        
        """
        if not hasattr(self, "_energy_max_database"):
            self._energy_max_database = float(max([m.energy for m in self.minima]))
        max_prob = float(self.minprob)
        energy_onset_width = 1.
        dE = self._energy_max_database + self.energy_offset - Emax
        f = dE / energy_onset_width
        if np.log(np.finfo('d').max) <= (-f):
            onset_prob = max_prob
        else:
            onset_prob = max_prob / (1. + np.exp(-f))
        return float(onset_prob)/self.nreplicas
    
    def get_starting_configurations(self, Emax):
        """this function overloads the function in NestedSampling"""
        # choose a replica randomly
        configs = self.get_starting_configurations_from_replicas()
        # replace each starting configuration with a one chosen
        # from the minima with probability onset_prob
        onset_prob = self.onset_prob_func(Emax)
        for i in range(len(configs)):
            if np.random.uniform(0,1) < onset_prob:
                x, energy = self.get_starting_configuration_minima(Emax)
                configs[i] = Replica(x, energy, from_random=False)
                if self.verbose:
                    print "sampling from minima, E minimum:", energy, "with probability:", onset_prob
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
    
