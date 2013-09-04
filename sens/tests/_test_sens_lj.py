import unittest
import numpy as np

from sens.bh_sampling import NestedSamplingSA
from sens.models._lj_tools import LJClusterSENS

from _test_ns_lj import TestNS_LJ

class TestSENS_LJ(TestNS_LJ):
    def setUp(self):
        self.setUp1()

    def setUp1(self, nproc=1):
        self.natoms = 13
        self.gmin = -44.3269
        self.system = LJClusterSENS(self.natoms, 2.5)
        self.nreplicas = 10
        self.stepsize = 0.01
        self.nproc = nproc
        
        self.database = self.system.create_database("lj13.db")
        self.minima = self.database.minima()
        
        self.mc_runner = self.system.get_mc_walker(mciter=100)

        self.ns = NestedSamplingSA(self.system, self.nreplicas, self.mc_runner,
                                   minima=self.minima, minprob=0.1,
                                   stepsize=0.1, nproc=nproc, verbose=True)
        
        self.Emax0 = self.ns.replicas[-1].energy
        
        self.niter = 100
        for i in xrange(self.niter):
            self.ns.one_iteration()
        self.Emax = self.ns.replicas[-1].energy
        self.Emin = self.ns.replicas[0].energy
    
    def test1(self):
        self.assert_(len(self.ns.replicas) == self.nreplicas)
        self.assert_(self.Emax < self.Emax0)
        self.assert_(self.Emin < self.Emax)
        self.assert_(self.Emin < 0)
        self.assert_(self.Emin >= self.gmin)
        self.assert_(self.ns.stepsize != self.stepsize)
        self.assertEqual(len(self.ns.max_energies), self.niter * self.nproc)

class testSENS_LJ_Par(TestSENS_LJ):
    def setUp(self):
        self.setUp1(nproc=3)
    
    
if __name__ == "__main__":
    unittest.main()  