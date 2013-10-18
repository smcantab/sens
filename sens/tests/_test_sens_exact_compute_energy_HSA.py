import unittest
import numpy as np
import itertools

from pele.mindist import TransformAtomicCluster
from pele.utils import rotations

from nested_sampling import Replica

from sens._HSA_sampler_cluster import HSASamplerCluster
from sens.models._lj_tools import LJClusterSENS
from _utils import build_database


class TestSENSExact_LJ(unittest.TestCase):
    def setUp(self):
        self.seed = np.random.randint(1000000)
#        self.seed = 549670 # failed self.assertGreater(self.ns.count_sampled_minima, 0)
        print "seed", self.seed
        np.random.seed(self.seed)
        self.setUp1()
    

    def setUp1(self, nproc=1):
        self.natoms = 31
        self.system = LJClusterSENS(self.natoms, 2.5)
        self.ndof = 3*self.natoms - 6

        self.nreplicas = 10
        self.stepsize = 0.01
        self.nproc = nproc
        
        
        try:
            self.database = self.system.create_database("/scratch/scratch2/js850/sens/testruns/lj31/lj31.sqlite", createdb=False)
        except IOError:
            self.database = build_database(self.system, 4)
            
        

        self.minima = list(self.database.minima())
        assert self.database.number_of_minima() > 1, "%d minima" %  self.database.number_of_minima()
        
        self.sampler = HSASamplerCluster(self.minima, self.ndof, copy_minima=True, 
                                         center_minima=True, 
                                         compare_structures=self.system.get_compare_exact(), 
                                         mindist=self.system.get_mindist(),
                                         minimizer=self.system.get_minimizer())

        self.Emax = min((m.energy for m in self.minima)) + 1.

        # for performing transformations like rotations, translations, permutations, etc
        self.transform = TransformAtomicCluster()
        
        # this is a list of random transformations that don't change the energy
        self.transformations = [self.rtrans, self.rrot, self.invert, self.rperm]

    def rtrans(self, x):
        vec3 = np.random.uniform(-1,1,3)
        self.transform.translate(x, vec3)
    
    def rrot(self, x):
        q = rotations.random_q()
        mx = rotations.q2mx(q)
        self.transform.rotate(x, mx)
    
    def invert(self, x):
        self.transform.invert(x)
    
    def rperm(self, x):
        perm = range(self.natoms)
        np.random.shuffle(perm)
        self.transform.permute(x, perm)
    
    def do_test(self, tform):
        # sample a configurations fr
        
        # sample a configuration from the HSA
        m, xsampled = self.sampler.sample_coords(self.Emax)
        
        # get the energy of that configuration in the HSA
        E_HSA = self.sampler.compute_energy(xsampled, m)
        
        # get the real energy of the sampled configuration
        pot = self.system.get_potential()
        Esampled = pot.getEnergy(xsampled)
        r = Replica(xsampled, Esampled)
        
        # transform xsampled in a way that preserves the real energy but not necessarily the HSA energy
        tform(xsampled)
        
        # assert the real energy has not changed
        Etrans = pot.getEnergy(xsampled)
        self.assertAlmostEqual(r.energy, Etrans, delta=1e-4)
        
        # use the nested sampling routine to compute the HSA energy of xsampled.
        # This should undo the transformations we applied and find the correct HSA energy 
        m_quench, E_HSA_computed = self.sampler.compute_energy_in_HSA(r.energy, r.x)
        
        # assert that self.ns._compute_energy_in_SA(r) was able to recover the correct HSA energy
        self.assertIsNotNone(E_HSA_computed)
        self.assertAlmostEqual(E_HSA, E_HSA_computed, delta=1e-4)
        self.assertAlmostEqual(m.energy, m_quench.energy, delta=1e-4)
        
    def fchain(self, flist):
        """return a function which applies all of the functions in flist to the input"""
        def function_chain(x):
            for f in reversed(flist):
                f(x)
        return function_chain
    
    def test1(self):
        """run do_test() on all combinations transformations"""
        for f in self.transformations:
            self.do_test(f)
    
    def test_2(self):
        """run do_test() on all combinations of length 2 of the transformations
        """
        for flist in itertools.product(self.transformations, repeat=2):
            tform = self.fchain(flist)
            self.do_test(tform)
    
    
    def test_10(self):
        """run do_test() on all combinations of length 10 of the transformations
        """
        maxiter = 500
        i = 0
        for flist in itertools.product(self.transformations, repeat=10):
            tform = self.fchain(flist)
            self.do_test(tform)
            i += 1
            if i >= maxiter: break
    
        
    
    
if __name__ == "__main__":
    unittest.main()  
