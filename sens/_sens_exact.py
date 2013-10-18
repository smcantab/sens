"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np
import time
import multiprocessing as mp
import collections
from nested_sampling import NestedSampling, Replica, Result

        
class _HSASwapper(object):
    """organize the swapping with the HSA
    """
    def __init__(self, hsa_sampler, potential, config_tests=None):
        self.hsa_sampler = hsa_sampler
        self.potential = potential
        self.config_tests = config_tests

    def attempt_swap(self, replica, Emax):
        res = Result()
        res.new_replica = None
        res.nfev = 0
        
        # sample a configuration from the harmonic superposition approximation
        m, xsampled = self.hsa_sampler.sample_coords(Emax)

        # if the configuration fails the config test then reject the swap
        if self.config_tests is not None:
            for test in self.config_tests:
                if not test(coords=xsampled):
                    res.reason_rejected = "config_test"
                    return res
        
        # if the energy returned by full energy function is too high, then reject the swap
        Esampled = self.potential.get_energy(xsampled)
        res.nfev += 1
        if Esampled >= Emax:
            res.reason_rejected = "Esampled too large"
            return res
        
        # compute the energy of the replica within the superposition approximation.
        # m_quench returns the minimum to which the replica is quenched to, 
        # thus identifying in which basin the replica is
        m_quench, E_HSA = self.hsa_sampler.compute_energy_in_HSA(replica.energy, replica.x)
        res.nfev += self.hsa_sampler._nfev_last_quench
        
        # reject if the returned energy is None.  This probably means that replica.x quenched to a minimum not in the database
        if E_HSA is None:
            res.reason_rejected = "E_HSA is None"
            return res

        # reject if the energy is too high
        if E_HSA >= Emax:
            # no swap done
            res.reason_rejected = "E_HSA too large"
            return res

        # return the results in a Results dictionary
        res.new_replica = Replica(xsampled, Esampled, from_random=False)
        res.minimum = m
        res.E_HSA = E_HSA
        res.m_quench = m_quench
        
        return res



class _HSASwapperParallel(mp.Process):
    def __init__(self, hsa_swapper, input_queue, output_queue):
        mp.Process.__init__(self)
        self.hsa_swapper = hsa_swapper
        self.input_queue = input_queue
        self.output_queue = output_queue
        
    def run(self):
        while True:
            val = self.input_queue.get()
            if val == "kill":
                break
            
            args, r_index = val
            result = self.hsa_swapper.attempt_swap(*args)
            self.output_queue.put((result, r_index))
        
    
class _SwapInfoAccumulator(object):
    def __init__(self):
        self._sampled_minima_counts = collections.Counter()
        self._effective_swap_counts_from = collections.Counter()
        self._effective_swap_counts_to = collections.Counter()
        self._reasons_rejected = collections.Counter()
        self.nfev = 0
    
    def add_swap(self, old_replica, result, Emax):
        self.nfev += result.nfev
        if result.new_replica is None:
            self._reasons_rejected[result.reason_rejected] += 1
        else:
            self._sampled_minima_counts[result.minimum.energy] += 1
            if result.m_quench.energy != result.minimum.energy:
                self._effective_swap_counts_from[result.m_quench.energy] += 1
                self._effective_swap_counts_to[result.minimum.energy] += 1
    
    def number_of_rejected_swaps(self):
        return sum(self._reasons_rejected.itervalues())
    
    def number_of_swaps(self):
        return sum(self._sampled_minima_counts.itervalues())
    
    def number_of_effective_swaps(self):
        return sum(self._effective_swap_counts_to.itervalues())
    
    def number_of_swap_attempts(self):
        return sum(self.number_of_rejected_swaps() + self.number_of_swaps())
    
    def minima_sampled_string(self):
        ostr = "minima sampled: "
        for energy, count in sorted(self._sampled_minima_counts.iteritems(), key=lambda mc: mc[0]):
            ostr += str(energy) + ": " + str(count) + ", "
        if len(self._sampled_minima_counts) > 0:
            # remove the comma and space before returning
            ostr = ostr[:-2]
        return ostr
    
    def effective_swap_to_minima_string(self):
        ostr = "effective swaps to minima: "
        for energy, count in sorted(self._effective_swap_counts_to.iteritems(), key=lambda mc: mc[0]):
            ostr += str(energy) + ": " + str(count) + ", "
        if len(self._effective_swap_counts_to) > 0:
            # remove the comma and space before returning
            ostr = ostr[:-2]
        return ostr
    
    def effective_swap_from_minima_string(self):
        ostr = "effective swaps from minima: "
        for energy, count in sorted(self._effective_swap_counts_from.iteritems(), key=lambda mc: mc[0]):
            ostr += str(energy) + ": " + str(count) + ", "
        if len(self._effective_swap_counts_from) > 0:
            # remove the comma and space before returning
            ostr = ostr[:-2]
        return ostr
    
    def reasons_rejected_string(self):
        ostr = "swaps rejected: "
        for reason, count in self._reasons_rejected.iteritems():
            ostr += str(reason) + ": " + str(count) + ", "
        if len(self._reasons_rejected) > 0:
            # remove the comma and space before returning
            ostr = ostr[:-2]
        return ostr
    
    def print_info(self):
        print self.minima_sampled_string()
        print self.effective_swap_to_minima_string()
        print self.effective_swap_from_minima_string()
        print self.reasons_rejected_string()
         
        
        
        


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
    def __init__(self, replicas, mc_runner,
                  hsa_sampler, potential, config_tests=None,
                  debug=True,
                  **kwargs):
        super(NestedSamplingSAExact, self).__init__(replicas, mc_runner, **kwargs)
        self.debug = debug
        self.config_tests = config_tests
        if self.debug and self.config_tests is None:
            print "warning: using no config tests"
        
        self.hsa_sampler = hsa_sampler
        self.hsa_swapper = _HSASwapper(self.hsa_sampler, potential, 
                                      config_tests=self.config_tests)
        
        self._times = Result(mc=0., sampling=0.)
        self._times.at_start = time.time()
        self._swap_info = _SwapInfoAccumulator()

        # Determine the minimum energy for swapping.
        # There is no point in swapping if there is only one minimum.
        # So stop swapping when Emax is less than the energy of the second lowest minimum
        energies = [m.energy for m in self.hsa_sampler.minima]
        energies.sort()
        self._swap_Emax_min = energies[1]
        if self.verbose:
            print "will stop swapping when Emax is less than", self._swap_Emax_min
        
        
        self._set_up_parallelization_swap()
        

    def number_swaps_accepted(self):
        return self._swap_info.number_of_swaps()

    def _set_up_parallelization_swap(self):
        """set up the parallel workers for doing swap attempts"""
        if self.nproc > 1:
            self._swap_get_queue = mp.Queue()
            self._swap_put_queue = mp.Queue()
            self._swap_workers = []
            for i in xrange(self.nproc):
                worker = _HSASwapperParallel(self.hsa_swapper, self._swap_put_queue, self._swap_get_queue)
                worker.daemon = True
                worker.start()
                self._swap_workers.append(worker)
    
    def finish(self):
        NestedSampling.finish(self)
        if self.nproc > 1:
            # terminate the swap workers
            for worker in self._swap_workers:
                self._swap_put_queue.put("kill")
            for worker in self._swap_workers:
                worker.join()
                worker.terminate()
                worker.join()

    
    def _accumulate_swap_info(self, old_replica, result, Emax):
        self._swap_info.add_swap(old_replica, result, Emax)
#        if result.new_replica is not None:
#            if self.verbose:
#                print "accepting swap %d: Eold %g Enew %g Eold_SA %g Emax %g m.energy %g" % (
#                            self.iter_number, old_replica.energy, result.new_replica.energy, 
#                            result.E_HSA, Emax, result.minimum.energy)


    def _attempt_swaps_parallel(self, replicas, Emax):
        """do all the swap attempts in parallel"""
        assert self.nproc > 1
        assert len(replicas) == self.nproc
        
        # put all the replicas into the queue
        for index, r in enumerate(replicas):
            self._swap_put_queue.put(((r, Emax), index))
        
        # retrieve the results from the queue
        for i in xrange(len(replicas)):
            res, index = self._swap_get_queue.get()
            self._accumulate_swap_info(replicas[index], res, Emax)
            if res.new_replica is not None:
#                self.count_sampled_minima += 1
                replicas[index] = res.new_replica
        
        assert self._swap_get_queue.empty() 

    def _attempt_swaps_serial(self, replicas, Emax):
        assert len(replicas) == 1
        """do all the swap attempts in serial""" 
        for i in xrange(len(replicas)):
            r = replicas[i]
            res = self.hsa_swapper.attempt_swap(r, Emax)
            self._accumulate_swap_info(replicas[i], res, Emax)
            if res.new_replica is not None:
#                self.count_sampled_minima += 1
                replicas[i] = res.new_replica

    def _attempt_swaps(self, replicas, Emax):
        """attempt to swap the replicas"""
        if Emax < self._swap_Emax_min:
            return
        if self.nproc > 1:
            return self._attempt_swaps_parallel(replicas, Emax)
        else:
            return self._attempt_swaps_serial(replicas, Emax)

    def do_monte_carlo_chain(self, replicas, Emax):
        """attempt to swap the replicas with one draw from the HSA before doing the monte carlo walk
        """
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
            print "times: mc walk", self._times.mc, "swapping", self._times.sampling, "tot", t2 - self._times.at_start, "this iter", t2-t1, t1-t0
            print "function evaluations: mc walk", self._mc_niter, "swapping", self._swap_info.nfev
            self._swap_info.print_info() 
        return replicas
        
        
