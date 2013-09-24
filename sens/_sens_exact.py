"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np
import time
import multiprocessing as mp

from nested_sampling import NestedSampling, Replica, Result

from sens._HSA_sampler_cluster import HSASamplerCluster


        
        


            


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
                  copy_minima=True, config_tests=None, minimizer=None, center_minima=False,
                  debug=True,
                  **kwargs):
        super(NestedSamplingSAExact, self).__init__(system, nreplicas, mc_runner, **kwargs)
        self.debug = debug
#        if copy_minima or center_minima:
#            self.minima = [_UnboundMinimum(m) for m in minima]
#            if center_minima:
#                # js850: this is a quick hack.  this should be done more elegantly
#                for m in self.minima:
#                    x = m.coords.reshape([-1,3])
#                    com = x.sum(0) / x.shape[0]
#                    x = x - com[np.newaxis,:]
#                    x = x.reshape(-1)
#                    m.coords[:] = x[:] 
#        else:
#            self.minima = minima
#        if self.verbose:
#            self._check_minima()
        
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
        
#        self.minima_searcher = _MinimaSearcher(self.minima, energy_accuracy, self.compare_structures)
        self.hsa_sampler = HSASamplerCluster(minima, self.system.k, copy_minima=copy_minima, 
                                             center_minima=center_minima, energy_accuracy=energy_accuracy, 
                                             compare_structures=self.compare_structures, 
                                             mindist=self.mindist, 
                                             minimizer=self.minimizer, 
                                             debug=self.debug)
        
        self.count_sampled_minima = 0
        self._times = Result(mc=0., sampling=0.)
        self._times.at_start = time.time()
        


    def _attempt_swap(self, replica, Emax):
        # sample a configuration from the harmonic superposition approximation
        m, xsampled = self.hsa_sampler.sample_coords(Emax)

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
        E_SA = self.hsa_sampler.compute_energy_in_HSA(replica.energy, replica.x)
        
        # reject if the energy is too high
        if E_SA is None or E_SA >= Emax:
            # no swap done
            return None

        if self.verbose:
            print "accepting swap %d: Eold %g Enew %g Eold_SA %g Emax %g m.energy %g" % (self.iter_number, replica.energy, Esampled, E_SA, Emax, m.energy)
        self.count_sampled_minima += 1
        
        return Replica(xsampled, Esampled, from_random=False)

#    def _attempt_swaps_parallel(self, replicas, Emax):
#        pool = mp.Pool(self.nproc)
#        results = pool.map(lambda r: self._attempt_swap(r, Emax), replicas)
#        for i in xrange(len(replicas)):
#            newr = results[i]
#            if newr is not None:
#                replicas[i] = newr

    def _attempt_swaps(self, replicas, Emax):
        for i in xrange(len(replicas)):
            r = replicas[i]
            rnew = self._attempt_swap(r, Emax)
            if rnew is not None:
                replicas[i] = rnew

    def do_monte_carlo_chain(self, replicas, Emax):
#        replicas = super(NestedSamplingSAExact, self).do_monte_carlo_chain(replicas, Emax)
        
        # try to swap this configuration with one sampled from the HSA
        t0 = time.time()
        self._attempt_swaps(replicas, Emax)
        
        # do a monte carlo walk
        t1 = time.time()
        replicas = super(NestedSamplingSAExact, self).do_monte_carlo_chain(replicas, Emax)
        
        t2 = time.time()
        self._times.mc += t2 - t1
        self._times.sampling += t1 - t0
        if self.verbose and self.iprint > 0 and self.iter_number % self.iprint == 0:
            print "time in mc walk", self._times.mc, "sampling", self._times.sampling, "tot", t2 - self._times.at_start, "this iter", t2-t1, t1-t0 
        return replicas
        
        
