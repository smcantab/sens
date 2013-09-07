import numpy as np
import sys

#from pele.utils.xyz import write_xyz
from pele.accept_tests import SphericalContainer
from pele.systems import LJCluster
from pele.mindist import PointGroupOrderCluster
from pele.optimize import Result

from nested_sampling.utils.rotations import vector_random_uniform_hypersphere

try:
    from sens.src.runmc import mc_cython as lj_mc_cython
except ImportError:
    print "warning, can't import compiled mc"


class SphericalContainerWraper(SphericalContainer):
    """return false if an atom is outside the spherical container"""
    def __call__(self, energy=None, coords=None):
        assert coords is not None
        return self.accept(coords)

class LJClusterSENS(LJCluster):
    """same as LJCluster, but attach some additional information"""
    def __init__(self, natoms, bradius):
        super(LJClusterSENS, self).__init__(natoms)
        self.nzero_modes = 6
        self.k = 3 * natoms - self.nzero_modes
        self.radius = bradius # 2.5 specific to lj31, 3 specific to lj38
        
        self.potential = self.get_potential()
    
    def get_energy(self, x):
        return self.potential.getEnergy(x)
    
#    def get_metric_tensor(self):
#        return None
    
#    def get_pgorder(self):
#        return PointGroupOrderCluster(self.get_compare_exact())
    
    def get_config_tests(self):
        return [SphericalContainerWraper(self.radius, nocenter=True)]
    
    def get_random_configuration(self):
        """make sure they're all inside the radius"""
        coords = np.zeros([self.natoms,3])
        for i in range(self.natoms):
            coords[i,:] = vector_random_uniform_hypersphere(3) * self.radius
            assert(np.linalg.norm(coords[i,:]) <= self.radius)

        # test to make sure the configuration is good
        tests = self.get_config_tests()
        for test in tests:
            assert test.accept(coords.flatten())
        return coords.flatten()

    def center_coords(self, x):
        x = x.reshape(-1,3)
        natoms = x.shape[0] 
        com = np.sum(x, 0) / natoms
        x -= com[np.newaxis, :]
    
    def get_mc_walker(self, mciter):
        return LJMonteCarloCompiled(self.radius, mciter=mciter)


class LJMonteCarloCompiled(object):
    def __init__(self, radius, mciter=100):
        self.radius = radius
        self.mciter = mciter
    
    def __call__(self, x0, stepsize, Emax, energy, seed=None):
        if seed is None:
            seed = np.random.randint(0, sys.maxint)
        x, energy, naccept = lj_mc_cython(x0, self.mciter, stepsize, Emax, self.radius, seed)
#        print ret
        res = Result()
        res.x0 = x0
        res.x = x
        res.nsteps = self.mciter
        res.naccept = naccept
        res.energy = energy
        return res
