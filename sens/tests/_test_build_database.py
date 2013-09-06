import unittest
import numpy as np

from pele.systems import LJCluster
from pele.thermodynamics import get_thermodynamic_information

from sens import get_all_normalmodes, NormalModes

class TestBuildDatabase(unittest.TestCase):
    def setUp(self):
        self.natoms = 6
        self.system = LJCluster(self.natoms)
        
        # create a database
        self.database = self.system.create_database()
        
        # add some minima to the database
        bh = self.system.get_basinhopping(self.database, outstream=None)
        while self.database.number_of_minima() < 2:
            bh.run(1)
        
        get_thermodynamic_information(self.system, self.database)
        

    def test(self):
        self.assertGreaterEqual(self.database.number_of_minima(), 2)
        for m in self.database.minima():
            self.assertIsNotNone(m.energy)
            self.assertIsNotNone(m.pgorder)
            self.assertIsNotNone(m.fvib)
        # test that the lowest minimum has the right vibrational free energy
        self.assertAlmostEqual(self.database.minima()[0].fvib, 56.085, 2)


class TestGetNormalModes(TestBuildDatabase):
    def setUp(self):
        super(TestGetNormalModes, self).setUp()
        
        get_all_normalmodes(self.system, self.database)
    
    def test(self):
        nzero = 6
        for m in self.database.minima():
            nm = m.normal_modes
            n = nm.freqs.size
            self.assertEqual((n, n), nm.vectors.shape)
            self.assertEqual(m, nm.minimum)
            for i in range(nzero):
                self.assertAlmostEqual(nm.freqs[i], 0., 3)
            self.assertGreater(nm.freqs[nzero], .01)
            
            v1 = nm.vectors[:,10]
            v2 = nm.vectors[:,11]
            dot = np.dot(v1, v2)
            self.assertAlmostEqual(dot, 0., 3)
            
            self.assertAlmostEqual(np.linalg.norm(v1), 1., 3)
            self.assertAlmostEqual(np.linalg.norm(v2), 1., 3)
            
    def test2(self):
        """compute the energies directly and with the SA and see that they agree"""
        nzero = 6
        pot = self.system.get_potential()
        for m in self.database.minima():
            nm = m.normal_modes
            
            for mode in xrange(nzero, nm.freqs.size):
                v = nm.vectors[:,mode].copy()
                freq = nm.freqs[mode]
                
                dist = 0.001
                x = m.coords + dist * v
                
                dE_computed = pot.getEnergy(x) - m.energy
                dE_SA = 0.5 * dist**2 * freq
                
#                print dE_computed, dE_SA
                
                self.assertAlmostEqual(dE_computed, dE_SA, 5)
            
            
             
            
        

if __name__ == "__main__":
    unittest.main()  
        