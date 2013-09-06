"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np
import bisect

from pele.storage.database import Minimum
from nested_sampling import NestedSampling, Replica

from sens import SASampler

class _UnboundMinimum(object):
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords
        

class _MinimaSearcher(object):
    def __init__(self, minima, energy_accuracy, compare_minima=None):
        self.minima = minima
        self.energy_accuracy = energy_accuracy
        self.compare_minima = compare_minima
        
        self.minima_dict = dict([(m.energy, m) for m in self.minima])
        self.energies = [m.energy for m in self.minima]
        self.energies.sort()
        
    def get_minima(self, energy, coords):
        """return the minimum that matches energy and coords.  return None if there is no match"""
        # first get a list of minima energies that might be the same minima
        Emax = energy + self.energy_accuracy
        Emin = energy - self.energy_accuracy
        imax = bisect.bisect_right(self.energies, Emax)
        imin = bisect.bisect_left(self.energies, Emin)
        
#        print "found", imax - imin, "cadidates"
        m = None
        mtest = _UnboundMinimum(energy, coords)
        for e in self.energies[imin:imax]:
            assert Emin <= e <= Emax
            m = self.minima_dict[e]

            # compare the coordinates more carefully
            if self.compare_minima is not None:            
                if self.compare_minima(mtest, m):
                    break
        
        return m
            

class ConfigTestError(StandardError):
    pass

class NestedSamplingSAExact(NestedSampling):
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
    def __init__(self, system, nreplicas, mc_runner, 
                  minima, energy_accuracy, compare_minima=None,
                  **kwargs):
        super(NestedSamplingSAExact, self).__init__(system, nreplicas, mc_runner, **kwargs)
        self.minima = minima
        self.minima_searcher = _MinimaSearcher(minima, energy_accuracy, compare_minima)
        
        self.sa_sampler = SASampler(self.minima, self.system.k)
        

        self.pot = self.system.get_potential()
        self.minimizer = self.system.get_minimizer()
        
        self.count_sampled_minima = 0



    def _compute_energy_in_SA(self, replica):
        # quench to nearest minimum
        qresult = self.minimizer(replica.x)
        
        # check if that minimum is in the database.  reject if not
        m = self.minima_searcher.get_minima(qresult.energy, qresult.coords)
        if m is None:
            return None
        
        
        
        pass

    def _attempt_swap(self, replica, Emax):
        m, xsampled = self.sa_sampler.sample_coords(Emax)
        Esampled = self.pot.get_energy(xsampled)
        
        if Esampled >= Emax:
            # no swap done
            return replica
        
        E_SA = self._compute_energy_in_SA(replica)
        
        if E_SA >= Emax:
            # no swap done
            return replica
        
        return Replica(xsampled, Esampled, from_random=True)
        

    def do_monte_carlo_chain(self, replicas, Emax):
        replicas = super(NestedSamplingSAExact, self).do_monte_carlo_chain(replicas, Emax)
        
        replicas = super(NestedSamplingSAExact, self).do_monte_carlo_chain(replicas, Emax)
        
        return replicas
        
        
