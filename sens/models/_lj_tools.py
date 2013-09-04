import numpy as np
import sys

#from pele.utils.xyz import write_xyz
from pele.accept_tests import SphericalContainer
from pele.systems import LJCluster
from pele.mindist import PointGroupOrderCluster
from pele.optimize import Result

from nested_sampling.utils.rotations import vector_random_uniform_hypersphere


class SphericalContainerWraper(SphericalContainer):
    """return false if an atom is outside the spherical container"""
    def __call__(self, energy=None, coords=None):
        return self.accept(coords)

class LJClusterNew(LJCluster):
    """same as LJCluster, but attach some additional information"""
    def __init__(self, natoms, bradius):
        super(LJClusterNew, self).__init__(natoms)
        self.nzero_modes = 6
        self.k = 3 * natoms - self.nzero_modes
        self.radius = bradius # 2.5 specific to lj31, 3 specific to lj38
    
#    def get_metric_tensor(self):
#        return None
    
    def get_pgorder(self):
        return PointGroupOrderCluster(self.get_compare_exact())
    
    def get_config_tests(self):
        return [SphericalContainerWraper(self.radius, nocenter=True)]
    
    def get_random_configuration(self):
        """make sure they're all inside the radius"""
        from pele.accept_tests import SphericalContainer
        test = self.get_config_tests()[0]
        coords = np.zeros([self.natoms,3])
        for i in range(self.natoms):
            coords[i,:] = vector_random_uniform_hypersphere(3) * self.radius
            assert(np.linalg.norm(coords[i,:]) <= self.radius)
        assert(test.accept(coords.flatten()))
        return coords.flatten()

    def center_coords(self, x):
        x = x.reshape(-1,3)
        natoms = x.shape[0] 
        com = np.sum(x, 0) / natoms
        x -= com[np.newaxis, :]


class MonteCarloCompiled(object):
    def __init__(self, radius):
        self.radius = radius
    
    def __call__(self, x0, mciter, stepsize, Emax, energy, seed=None):
        if seed is None:
            seed = np.random.randint(0, sys.maxint)
        x, energy, naccept = mc_cython(x0, mciter, stepsize, Emax, self.radius, seed)
#        print ret
        res = Result()
        res.x0 = x0
        res.x = x
        res.nsteps = mciter
        res.naccept = naccept
        res.energy = energy
        return res
