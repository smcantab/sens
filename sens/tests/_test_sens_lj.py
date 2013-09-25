import unittest
import numpy as np

from nested_sampling import Replica
from sens import NestedSamplingSA

import _test_ns_lj
from _utils import build_database, create_replicas

class TestSENS_LJ(_test_ns_lj.TestNS_LJ):
    def setUp(self):
        self.setUp1()

    def setUp1(self, nproc=1):
        self.set_up_system()
        self.nreplicas = 10
        self.stepsize = 0.01
        self.nproc = nproc
        
        try:
            self.database = self.system.create_database("lj13.db", createdb=False)
        except IOError:
            self.database = build_database(self.system, 20, dbfname="lj13.db")
            
        self.minima = self.database.minima()
        assert self.database.number_of_minima() > 1, "%d minima" %  self.database.number_of_minima()
        
        self.mc_runner = self.system.get_mc_walker(mciter=100)

        replicas = create_replicas(self.system, self.nreplicas)
        self.ns = NestedSamplingSA(replicas, self.mc_runner, self.minima, self.system.k,
                                   config_tests=self.system.get_config_tests(),
                                   minprob=0.9, energy_offset=100.,
                                   stepsize=0.1, nproc=nproc, verbose=False,
                                   center_minima=True)
        
        self.Emax0 = self.ns.replicas[-1].energy
        
        self.niter = 100
        for i in xrange(self.niter):
            self.ns.one_iteration()
        self.Emax = self.ns.replicas[-1].energy
        self.Emin = self.ns.replicas[0].energy
    
    def test2(self):
        self.assertGreater(self.ns.count_sampled_minima, 0)
    

class testSENS_LJ_Par(TestSENS_LJ):
    def setUp(self):
        self.setUp1(nproc=3)
    
    
if __name__ == "__main__":
    unittest.main()  
