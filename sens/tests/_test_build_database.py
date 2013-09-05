import unittest

from pele.systems import LJCluster
from pele.thermodynamics import get_thermodynamic_information

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

    
if __name__ == "__main__":
    unittest.main()  
        