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
                  minima, energy_accuracy, compare_minima=None, mindist=None,
                  **kwargs):
        super(NestedSamplingSAExact, self).__init__(system, nreplicas, mc_runner, **kwargs)
        self.minima = minima
        if self.verbose:
            self._check_minima()
        self.minima_searcher = _MinimaSearcher(minima, energy_accuracy, compare_minima)
        self.mindist = mindist
        
        self.sa_sampler = SASampler(self.minima, self.system.k)
        

        self.minimizer = self.system.get_minimizer()
        
        self.count_sampled_minima = 0
        


    def _check_minima(self):
        for m in self.minima:
            assert m.energy is not None
            assert m.fvib is not None
            assert m.pgorder is not None
            assert m.normal_modes is not None

    def _compute_energy_in_SA(self, replica):
        # quench to nearest minimum
        qresult = self.minimizer(replica.x)
        
        # check if that minimum is in the database.  reject if not
        m = self.minima_searcher.get_minima(qresult.energy, qresult.coords)
        if m is None:
            return None
        
        # put the two structures in best alignment.
        # e.g. account for trivial tranlational and rotational degrees of freedom
        x, x0 = replica.x.copy(), qresult.coords.copy()
        if self.mindist is not None:
            dist, x, x0 = self.mindist(x, x0)
        
        energy = self.sa_sampler.compute_energy(x, x0, m)
        
        return energy

    def _attempt_swap(self, replica, Emax):
        # sample a configuration from the harmonic superposition approximation
        m, xsampled = self.sa_sampler.sample_coords(Emax)
        
        # if the energy returned by full energy function is too high, then reject the swap
        Esampled = self.system.get_energy(xsampled)
        if Esampled >= Emax:
            return replica
        
        # compute the energy of the replica within the superposition approximation.
        E_SA = self._compute_energy_in_SA(replica)
        
        # reject if the energy is too high
        if E_SA is None or E_SA >= Emax:
            # no swap done
            return replica

        print "accepting swap"
        self.count_sampled_minima += 1
        
        return Replica(xsampled, Esampled, from_random=True)
        

    def do_monte_carlo_chain(self, replicas, Emax):
        replicas = super(NestedSamplingSAExact, self).do_monte_carlo_chain(replicas, Emax)
        
        print "attempting swaps"
        for i, r in enumerate(replicas):
            rnew = self._attempt_swap(r, Emax)
            if rnew is not None:
                replicas[i] = rnew
        
        replicas = super(NestedSamplingSAExact, self).do_monte_carlo_chain(replicas, Emax)
        
        return replicas
        
        
