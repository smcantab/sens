"""
this test will use the harmonic superposition approximation generate random configurations for the nested sampling algorithm.  
This should be exactly equal to the directly computable HSA results 
"""
import unittest
import numpy as np

from pele.optimize import Result
from pele.thermodynamics import minima_to_cv

from nested_sampling import NestedSampling, run_nested_sampling, Replica, compute_heat_capacity

from sens.tests._utils import build_database
from sens.models import LJClusterSENS
from sens import SASampler

    

class SamplerSystem(object):
    def __init__(self, natoms=6):
        self.system = LJClusterSENS(natoms, 5.)
        self.ndof = natoms * 3 - 6
        self.natoms = natoms
        
        self.database = build_database(self.system, nminima=2)
        
        self.sampler = SASampler(self.database.minima(), self.ndof)
    
    def sample_configuration(self, Emax):
        m, x = self.sampler.sample_coords(Emax)
        E = self.sampler.compute_energy(x, m)
        return E, x
    
    def __call__(self, *args, **kwargs):
        Emax = args[2]
#        print "walking with Emax", Emax
        E, x = self.sample_configuration(Emax)
        
        res = Result(x=x, energy=E, naccept=5, nsteps=10, Emax=Emax)
        return res
        
    
class TestSamplerNS(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        self.npar = 1
        self.system = SamplerSystem()
        
        
        self.nreplicas = 200
        replicas = [Replica(x, energy) for energy, x in 
                    [self.system.sample_configuration(10.) for i in range(self.nreplicas)]]
        
        mc_walker = self.system 
        ns = NestedSampling(self.system, nreplicas=self.nreplicas, 
                            mc_walker=mc_walker, nproc=1, verbose=False, 
                            iprint=100, replicas=replicas)
        run_nested_sampling(ns, label="test")
        
        self.energies = np.array(ns.max_energies)
        
        self.T = np.linspace(0.1, 1.0, 10)
        # compute the heat capacity and internal energy from the results of the NS run 
        self.cvNS = compute_heat_capacity(self.energies, self.nreplicas, npar=self.npar, 
                                       ndof=self.system.ndof, 
                                       Tmin=min(self.T), Tmax=max(self.T), nT=len(self.T), 
                                       live_replicas=False)
        
        # compute the same directly from the HSA
        self.cvHSA = minima_to_cv(self.system.database.minima(), kT=self.T, k=self.system.ndof)
        

    def test1(self):
        for i in range(len(self.T)):
            self.assertAlmostEqual(self.cvNS.U[i], self.cvHSA.U[i], delta=.5)
            self.assertAlmostEqual(self.cvNS.U2[i], self.cvHSA.U2[i], delta=2.5)
            self.assertAlmostEqual(self.cvNS.Cv[i], self.cvHSA.Cv[i], delta=2.5)
        

if __name__ == "__main__":
    unittest.main()