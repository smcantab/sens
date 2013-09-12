import unittest
import numpy as np

from sens import NestedSamplingSA

import _test_ns_lj
from _utils import build_database

class TestSENS_LJ_Long(_test_ns_lj.TestNS_LJ):
    def setUp(self):
        self.setUp1()

    def setUp1(self, nproc=4):
        self.set_up_system()
        self.nreplicas = 50
        self.stepsize = 0.01
        self.nproc = nproc

        try:
            self.database = self.system.create_database("lj13.db", createdb=False)
        except IOError:
            self.database = build_database(self.system, 20, dbfname="lj13.db")
        self.minima = self.database.minima()

        self.mc_runner = self.system.get_mc_walker(mciter=10000)

        self.ns = NestedSamplingSA(self.system, self.nreplicas, self.mc_runner,
                                   minima=self.minima, minprob=0.1,
                                   stepsize=0.1, nproc=nproc, verbose=True,
                                   iprint=100)
        
        self.Emax0 = self.ns.replicas[-1].energy
        
        max_iter = 10000
        self.Etol = .01
        for i in xrange(max_iter):
            self.ns.one_iteration()
            deltaE = self.ns.replicas[-1].energy - self.ns.replicas[0].energy
            if  deltaE < self.Etol:
                break
        self.niter = i + 1
        self.Emax = self.ns.replicas[-1].energy
        self.Emin = self.ns.replicas[0].energy
    

    def test1(self):
        super(TestSENS_LJ_Long, self).test1()
        self.assertTrue(self.Emin < self.gmin + 1.,
                        "Nested sampling did not get to the bottom of the landscape: %g != %g" % (self.gmin, self.Emin))
        self.assertGreater(self.ns.count_sampled_minima, 0)

    
if __name__ == "__main__":
    unittest.main()  
