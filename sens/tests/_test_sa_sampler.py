import unittest
import numpy as np

from pele.systems import LJCluster
from pele.thermodynamics import get_thermodynamic_information

from sens import SASampler
from sens import get_all_normalmodes

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
        get_all_normalmodes(self.system, self.database)
        
        self.ndof = self.natoms * 3 - 6
        self.sampler = SASampler(self.database.minima(), self.ndof)

#    def test(self):
#        for m in self.database.minima():
#            print m.energy
#        print self.sampler.compute_weights(-11.7)
#        for i in range(10):
#            m = self.sampler.sample_minimum(-11)
#            print m._id

    def test1(self):
        Emax = -11.7
        minima, weights = self.sampler.compute_weights(0.)
        minima, weights = self.sampler.compute_weights(Emax)
        self.assertAlmostEqual(weights[0], 0.81256984, 3)
        self.assertAlmostEqual(weights[1], 1., 7)

    def test2(self):
        Emax = -11.7
        minima, weights = self.sampler.compute_weights(Emax)
        minima, weights = self.sampler.compute_weights(Emax)
        self.assertAlmostEqual(weights[0], 0.81256984, 3)
        self.assertAlmostEqual(weights[1], 1., 7)
    
    def test3(self):
        m = self.sampler.sample_minimum(-11.7)
        self.assertIn(m, self.database.minima())
    
    def test4(self):
        pot = self.system.get_potential()
        for m in self.database.minima():
            x = m.coords.copy()
            x += np.random.uniform(-1e-3, 1e-3, x.shape)
            ehsa = self.sampler.compute_energy(x, m.coords, m)
            ecalc = pot.getEnergy(x)
            ecompare = (ehsa - ecalc) / (ecalc - m.energy) 
            print ehsa - m.energy, ecalc - m.energy, m.energy, ecompare
            self.assertAlmostEqual(ecompare, 0., 1)
            

    
#    def test4(self):
#        Emax = -11.7
#        m = self.sampler.sample_minimum(Emax)
#        coords = self.sampler.sample_coords_from_basin(m, Emax)
#        self.assertEqual(coords.size, 3*self.natoms)
        

if __name__ == "__main__":
    unittest.main()  
