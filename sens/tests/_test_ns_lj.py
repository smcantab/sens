import unittest
import numpy as np

from nested_sampling import NestedSampling
from nested_sampling.src.cv_trapezoidal import compute_cv_c

from sens.models._lj_tools import LJClusterSENS

class TestNS_LJ(unittest.TestCase):
    def setUp(self):
        self.setUp1()

    def compute_cv(self, Tmin=.01, Tmax=3., nT=100):
        T = np.arange(Tmin, Tmax, (Tmax - Tmin) / nT)
        Cv = compute_cv_c(np.array(self.ns.max_energies), 
                          float(self.nproc), float(self.ndof), T.min(), T.max(), T.size, float(self.ndof), live=False)
        return T, Cv


    def set_up_system(self):
        self.natoms = 13
        self.gmin = -44.3269
        self.system = LJClusterSENS(self.natoms, 2.5)
        self.ndof = 3 * self.natoms - 6

    def setUp1(self, nproc=1):
        self.set_up_system()
        self.nreplicas = 10
        self.stepsize = 0.01
        self.nproc = nproc
        
        self.mc_runner = self.system.get_mc_walker(mciter=100)

        self.ns = NestedSampling(self.system, self.nreplicas, self.mc_runner, 
                                 stepsize=0.1, nproc=nproc, verbose=False)
        
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
        self.assertEqual(self.ns.failed_mc_walks, 0)

    def run_ns(self, max_iter=100, Etol=1e-4):
#        max_iter = 10000
#        self.Etol = .01
        self.Etol = Etol
        for i in xrange(int(max_iter)):
            self.ns.one_iteration()
            deltaE = self.ns.replicas[-1].energy - self.ns.replicas[0].energy
            if  deltaE < self.Etol:
                break
        self.niter = i + 1
        self.Emax = self.ns.replicas[-1].energy
        self.Emin = self.ns.replicas[0].energy

#    def testcv(self):
#        T, cv = self.compute_cv()
#        import matplotlib.pyplot as plt
##        plt.plot(cv)
#        print cv.shape
#        plt.plot(T, cv)
#        plt.show()


class testNSPar(TestNS_LJ):
    def setUp(self):
        self.setUp1(nproc=3)
    
    
if __name__ == "__main__":
    unittest.main()  