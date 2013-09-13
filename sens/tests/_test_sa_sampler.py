import unittest
import numpy as np

from pele.systems import LJCluster
from pele.thermodynamics import get_thermodynamic_information
import pele.utils.rotations as rotations


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

    def compute_weights(self, Emax, k):
        print [m.energy for m in self.database.minima()], Emax
        lweights = [ - np.log(m.pgorder) + 0.5 * k * np.log(Emax - m.energy) - 0.5 * m.fvib
                    for m in self.database.minima() if m.energy < Emax]
        lweights = np.array(lweights)
        if lweights.size <= 1: return lweights
        lwmax = lweights.max()
        lweights -= lwmax
        return np.exp(lweights)

    def test1(self):
        Emax = -11.7
        minima, weights = self.sampler.compute_weights(0.)
        minima, weights = self.sampler.compute_weights(Emax)
        self.assertAlmostEqual(weights[0], 0.81256984, 3)
        self.assertAlmostEqual(weights[1], 1., 7)
        
        new_weights = self.compute_weights(Emax, self.ndof)
        print weights
        print new_weights
        

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
        """check that the HSA energy is computed correctly"""
        pot = self.system.get_potential()
        for m in self.database.minima():
            x = m.coords.copy()
            x += np.random.uniform(-1e-3, 1e-3, x.shape)
            ehsa = self.sampler.compute_energy(x, m, x0=m.coords)
            ecalc = pot.getEnergy(x)
            ecompare = (ehsa - ecalc) / (ecalc - m.energy) 
            print ehsa - m.energy, ecalc - m.energy, m.energy, ecompare
            self.assertAlmostEqual(ecompare, 0., 1)
    
    def test_rotation(self):
        """assert that the HSA energy is *not* invariant under rotation"""
        pot = self.system.get_potential()
        aa = rotations.random_aa()
        rmat = rotations.aa2mx(aa)
        from pele.mindist import TransformAtomicCluster
        tform = TransformAtomicCluster(can_invert=True)
        for m in self.database.minima():
            x = m.coords.copy()
            # randomly move the atoms by a small amount
            x += np.random.uniform(-1e-3, 1e-3, x.shape)
            ehsa1 = self.sampler.compute_energy(x, m, x0=m.coords)
            ecalc = pot.getEnergy(x)
            # now rotate by a random matrix
            xnew = x.copy()
            tform.rotate(xnew, rmat)
            ehsa2 = self.sampler.compute_energy(xnew, m, x0=m.coords)
            ecalc2 = pot.getEnergy(xnew)
            self.assertAlmostEqual(ecalc, ecalc2, 5)
            self.assertNotAlmostEqual(ehsa1, ehsa2, 1)
    
#    def test_rotation_2(self):
#        """assert that the HSA energy *is* invariant under rotation *if* the initial coords are also rotated"""
#        pot = self.system.get_potential()
#        aa = rotations.random_aa()
#        rmat = rotations.aa2mx(aa)
#        from pele.mindist import TransformAtomicCluster
#        tform = TransformAtomicCluster(can_invert=True)
#        for m in self.database.minima():
#            x = m.coords.copy()
#            # randomly move the atoms by a small amount
#            x += np.random.uniform(-1e-3, 1e-3, x.shape)
#            ehsa1 = self.sampler.compute_energy(x, m, x0=m.coords)
#            ecalc = pot.getEnergy(x)
#            # now rotate by a random matrix
#            xnew = x.copy()
#            tform.rotate(xnew, rmat)
#            xmnew = m.coords.copy()
#            tform.rotate(xmnew, rmat)
#            ehsa2 = self.sampler.compute_energy(xnew, m, x0=xmnew)
#            ecalc2 = pot.getEnergy(xnew)
#            self.assertAlmostEqual(ecalc, ecalc2, 5)
#            self.assertAlmostEqual(ehsa1, ehsa2, 3)
    
    def test_permutation(self):
        """assert that the HSA energy is not invariant under permutation"""
        pot = self.system.get_potential()
        perm = range(self.natoms)
        np.random.shuffle(perm)
        from pele.mindist import TransformAtomicCluster
        tform = TransformAtomicCluster(can_invert=True)
        for m in self.database.minima():
            x = m.coords.copy()
            # randomly move the atoms by a small amount
            x += np.random.uniform(-1e-3, 1e-3, x.shape)
            ehsa1 = self.sampler.compute_energy(x, m, x0=m.coords)
            ecalc = pot.getEnergy(x)
            # now rotate by a random matrix
            xnew = x.copy()
            xnew = tform.permute(xnew, perm)
            ehsa2 = self.sampler.compute_energy(xnew, m, x0=m.coords)
            ecalc2 = pot.getEnergy(xnew)
            self.assertAlmostEqual(ecalc, ecalc2, 5)
            self.assertNotAlmostEqual(ehsa1, ehsa2, 1)
            
#    def test_permutation_2(self):
#        """assert that the HSA energy *is* invariant under permutation *if* the reference coords are also permuted"""
#        pot = self.system.get_potential()
#        perm = range(self.natoms)
#        np.random.shuffle(perm)
#        from pele.mindist import TransformAtomicCluster
#        tform = TransformAtomicCluster(can_invert=True)
#        for m in self.database.minima():
#            x = m.coords.copy()
#            # randomly move the atoms by a small amount
#            x += np.random.uniform(-1e-3, 1e-3, x.shape)
#            ehsa1 = self.sampler.compute_energy(x, m, x0=m.coords)
#            ecalc = pot.getEnergy(x)
#            # now rotate by a random matrix
#            xnew = x.copy()
#            xnew = tform.permute(xnew, perm)
#            xmnew = m.coords.copy()
#            xmnew = tform.permute(xmnew, perm)
#            ehsa2 = self.sampler.compute_energy(xnew, m, x0=xmnew)
#            ecalc2 = pot.getEnergy(xnew)
#            self.assertAlmostEqual(ecalc, ecalc2, 5)
#            self.assertAlmostEqual(ehsa1, ehsa2, 3)
#            
#    def test_rotation_permutation_2(self):
#        """assert that the HSA energy *is* invariant under permutation *if* the reference coords are also permuted"""
#        pot = self.system.get_potential()
#        aa = rotations.random_aa()
#        rmat = rotations.aa2mx(aa)
#
#        perm = range(self.natoms)
#        np.random.shuffle(perm)
#        from pele.mindist import TransformAtomicCluster
#        tform = TransformAtomicCluster(can_invert=True)
#        for m in self.database.minima():
#            x = m.coords.copy()
#            # randomly move the atoms by a small amount
#            x += np.random.uniform(-1e-3, 1e-3, x.shape)
#            ehsa1 = self.sampler.compute_energy(x, m, x0=m.coords)
#            ecalc = pot.getEnergy(x)
#            # now rotate by a random matrix
#            xnew = x.copy()
#            xnew = tform.permute(xnew, perm)
#            xmnew = m.coords.copy()
#            xmnew = tform.permute(xmnew, perm)
#            tform.rotate(xnew, rmat)
#            tform.rotate(xmnew, rmat)
#            ehsa2 = self.sampler.compute_energy(xnew, m, x0=xmnew)
#            ecalc2 = pot.getEnergy(xnew)
#            self.assertAlmostEqual(ecalc, ecalc2, 5)
#            self.assertAlmostEqual(ehsa1, ehsa2, 3)
            
#            ecompare = (ehsa - ecalc) / (ecalc - m.energy) 
#            print ehsa - m.energy, ecalc - m.energy, m.energy, ecompare
#            self.assertAlmostEqual(ecompare, 0., 1)

            

    
#    def test4(self):
#        Emax = -11.7
#        m = self.sampler.sample_minimum(Emax)
#        coords = self.sampler.sample_coords_from_basin(m, Emax)
#        self.assertEqual(coords.size, 3*self.natoms)
        

if __name__ == "__main__":
    unittest.main()  
