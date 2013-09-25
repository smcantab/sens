import unittest
import numpy as np

from nested_sampling import NestedSampling

import _test_ns_lj
import _utils

class TestNS_LJ_Long(_test_ns_lj.TestNS_LJ):
    def setUp(self):
        self.setUp1()

    def setUp1(self, nproc=4):
        self.set_up_system()
        self.nreplicas = 50
        self.stepsize = 0.01
        self.nproc = nproc
        
        self.mc_runner = self.system.get_mc_walker(mciter=10000)

        replicas = _utils.create_replicas(self.system, self.nreplicas)
        self.ns = NestedSampling(replicas, self.mc_runner, 
                                 stepsize=self.stepsize, nproc=nproc, 
                                 verbose=True, iprint=100)
        
        self.Emax0 = self.ns.replicas[-1].energy
        
        self.run_ns(max_iter=10000, Etol=.01)
#        max_iter = 10000
#        self.Etol = .01
#        for i in xrange(max_iter):
#            self.ns.one_iteration()
#            deltaE = self.ns.replicas[-1].energy - self.ns.replicas[0].energy
#            if  deltaE < self.Etol:
#                break
#        self.niter = i + 1
#        self.Emax = self.ns.replicas[-1].energy
#        self.Emin = self.ns.replicas[0].energy
    

    def test1(self):
        super(TestNS_LJ_Long, self).test1()
        self.assertTrue(self.Emin < self.gmin + 1.,
                        "Nested sampling did not get to the bottom of the landscape: %g != %g" % (self.gmin, self.Emin))

    
if __name__ == "__main__":
    unittest.main()  