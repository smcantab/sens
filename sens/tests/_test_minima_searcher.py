import unittest
import numpy as np

from sens._sens_exact import _MinimaSearcher
from sens.models._lj_tools import LJClusterSENS

class TestMinimaSearch(unittest.TestCase):
    def setUp(self):
        self.natoms = 6
        self.system = LJClusterSENS(self.natoms, 2.5)
        
        # create a database
        self.database = self.system.create_database()
        
        # add some minima to the database
        bh = self.system.get_basinhopping(self.database, outstream=None)
        while self.database.number_of_minima() < 2:
            bh.run(1)
        
        self.ndof = self.natoms * 3 - 6
    
        self.minima_searcher = _MinimaSearcher(self.database.minima(), energy_accuracy=1e-4, compare_minima=self.database.compareMinima)

    def test_exact(self):
        for m in self.database.minima():
            mret = self.minima_searcher.get_minima(m.energy, m.coords)
            self.assertEqual(m, mret)

    def test_random(self):
        for m in self.database.minima():
            mret = self.minima_searcher.get_minima(m.energy + np.random.uniform(-1e-4, 1e-4), m.coords)
            self.assertEqual(m, mret)
    def test_none(self):
        for m in self.database.minima():
            mret = self.minima_searcher.get_minima(m.energy + 10., m.coords)
            self.assertIsNone(mret)


if __name__ == "__main__":
    unittest.main()  
