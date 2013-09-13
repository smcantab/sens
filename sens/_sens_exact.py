"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np
import bisect
from collections import namedtuple

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
    compare_structures : callable, `compare_structures(coords1, coords2)`
        a function that returns true if two structures are the same and false otherwise.
        This is used only if the two minima have energies that are within 
        energy_accuracy of each other.
    """
    def __init__(self, minima, energy_accuracy, compare_structures=None):
        self.minima = minima
        self.energy_accuracy = energy_accuracy
        self.compare_structures = compare_structures
        
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
        RetVal = namedtuple("RetVal", ["minimum", "transformation"])
        retval = RetVal(minimum=None, transformation=None)
        mtest = _UnboundMinimumSmall(energy, coords)
        for e in self.energies[imin:imax]:
            assert Emin <= e <= Emax
            m = self.minima_dict[e]
            retval = RetVal(minimum=m, transformation=None)

            
            # compare the coordinates more carefully
            if self.compare_structures is not None:
                transform = self.compare_structures.find_transformation(m.coords, mtest.coords)
                if transform is not None:
                    retval._replace(transformation=transform)
                    break
        
        return retval
            

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
                  minima, energy_accuracy, compare_structures=None, mindist=None,
                  copy_minima=True, config_tests=None, minimizer=None,
                  debug=True,
                  **kwargs):
        super(NestedSamplingSAExact, self).__init__(system, nreplicas, mc_runner, **kwargs)
        self.debug = debug
        if copy_minima:
            self.minima = [_UnboundMinimum(m) for m in minima]
        else:
            self.minima = minima
        if self.verbose:
            self._check_minima()
        
        self.compare_structures = compare_structures
        if compare_structures is None:
            try:
                self.compare_structures = self.system.get_compare_exact()
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
        
        self.minima_searcher = _MinimaSearcher(self.minima, energy_accuracy, self.compare_structures)
        self.sa_sampler = SASampler(self.minima, self.system.k)

        
        self.count_sampled_minima = 0
        
    def _check_minima(self):
        for m in self.minima:
            assert m.energy is not None
            assert m.fvib is not None
            assert m.pgorder is not None
            assert m.normal_modes is not None

    def _compute_energy_in_SA(self, replica):
        """compute the energy of the coordinates in replica.coords in the harmonic superposition approximation
        
        This is where most of the difficulty is in the algorithm.  
        
        This is also where all the system dependence is
        """
        # quench to nearest minimum
        qresult = self.minimizer(replica.x)
        
        # check if that minimum is in the database.  reject if not
        m, transformation = self.minima_searcher.get_minima(qresult.energy, qresult.coords)
        if m is None:
            return None
        
        # put replica.coords into best alignment with the structure stored in m.coords
        # this involves accounting for symmetries of the Hamiltonian like translational, 
        # rotational and permutations symmetries.  You can use the coordinates in qresult.coords
        # to help find the best permutation.  Ultimately it must be aligned with the structure in m.coords
        # The hessian eigenvectors were computed with a given permutation
        # e.g. account for trivial translational and rotational degrees of freedom
        x = replica.x.copy()
        if transformation is not None:
            # transformation is the set of transformations that put qresult.coords into exact
            # alignment with m.coords.  If we apply these transformations to replica.x then
            # replica.x will be in good (although not perfect) alignment already with m.coords 
            x = self.compare_structures.apply_transformation(x, transformation)
            
        # Do a final round of optimization to further improve the alignment
        if self.mindist is not None:
            dist, x0, x = self.mindist(m.coords.copy(), x)
            if self.debug:
                diff = np.linalg.norm(x0 - m.coords)
                if diff > .01:
                    with open("error.xyz", "w") as fout:
                        from pele.utils.xyz import write_xyz
                        write_xyz(fout, x0)
                        write_xyz(fout, m.coords)
                        
                    raise Exception("warning, mindist appears to have changed x0.  the norm of the difference is %g" % diff)
                    
                assert (x0 == m.coords).all()
        
        energy = self.sa_sampler.compute_energy(x, m)
        
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
            print "accepting swap: Eold %g Enew %g Eold_SA %g Emax %g m.energy %g" % (replica.energy, Esampled, E_SA, Emax, m.energy)
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
        
        
