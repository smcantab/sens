"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np
import bisect

from nested_sampling import NestedSampling, Replica

from sens import SASampler

class _UnboundMinimumSmall(object):
    """represent a minimum object unbound from the database"""
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords

class _UnboundNormalModes(object):
    def __init__(self, nm):
        self.freqs = nm.freqs.copy()
        self.vectors = nm.vectors.copy()

class _UnboundMinimum(object):
    """represent a minimum object unbound from the database"""
    def __init__(self, m):
        self.energy = m.energy
        self.coords = m.coords
        self.fvib = m.fvib
        self.pgorder = m.pgorder
        self.normal_modes = _UnboundNormalModes(m.normal_modes)
        
        

class _MinimaSearcher(object):
    """class to manage searching efficiently for minima in a database
    
    Parameters
    ----------
    minima : list of Mimumum objects
    energy_accuracy : float
        energy tolerance for when to consider that two minima are the same
    compare_minima : callable, `compare_minima(min1, min2)`
        a function that returns true if two minima are the same and false otherwise.
        This is used only if the two minima have energies that are within 
        energy_accuracy of each other.
    """
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
        mtest = _UnboundMinimumSmall(energy, coords)
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
                  copy_minima=True, config_tests=None, minimizer=None,
                  **kwargs):
        super(NestedSamplingSAExact, self).__init__(system, nreplicas, mc_runner, **kwargs)
        if copy_minima:
            self.minima = [_UnboundMinimum(m) for m in minima]
        else:
            self.minima = minima
        if self.verbose:
            self._check_minima()
        
        self.compare_minima = compare_minima
        if compare_minima is None:
            try:
                self.compare_minima = self.system.get_compare_minima()
            except NotImplementedError or AttributeError:
                pass
        
        self.mindist = mindist
        if mindist is None:
            try:
                self.mindist = self.system.get_mindist()
            except NotImplementedError or AttributeError:
                pass

        self.config_tests = config_tests
        if config_tests is None:
            try:
                self.config_tests = self.system.get_config_tests()
            except NotImplementedError or AttributeError:
                pass

        self.minimizer = minimizer
        if self.minimizer is None:
            self.minimizer = self.system.get_minimizer()
        
        self.minima_searcher = _MinimaSearcher(minima, energy_accuracy, compare_minima)
        self.sa_sampler = SASampler(self.minima, self.system.k)

        
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

        # if the configuration fails the config test then reject the swap
#        print "attempting swap"
        if self.config_tests is not None:
            for test in self.config_tests:
                if not test(coords=xsampled):
                    return None
        
        # if the energy returned by full energy function is too high, then reject the swap
        Esampled = self.system.get_energy(xsampled)
        if Esampled >= Emax:
            return None
        
        # compute the energy of the replica within the superposition approximation.
        E_SA = self._compute_energy_in_SA(replica)
        
        # reject if the energy is too high
        if E_SA is None or E_SA >= Emax:
            # no swap done
            return None

        if self.verbose:
            print "accepting swap: Eold %g Enew %g Eold_SA %g Emax %g" % (replica.energy, Esampled, E_SA, Emax)
        self.count_sampled_minima += 1
        
        return Replica(xsampled, Esampled, from_random=False)
        

    def do_monte_carlo_chain(self, replicas, Emax):
#        replicas = super(NestedSamplingSAExact, self).do_monte_carlo_chain(replicas, Emax)
        
        # try to swap this configuration with one sampled from the HSA
        for i in xrange(len(replicas)):
            r = replicas[i]
            rnew = self._attempt_swap(r, Emax)
            if rnew is not None:
                replicas[i] = rnew
        
        # do a monte carlo walk
        replicas = super(NestedSamplingSAExact, self).do_monte_carlo_chain(replicas, Emax)
        
        return replicas
        
        
