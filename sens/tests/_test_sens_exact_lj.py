import unittest
import numpy as np

from pele.thermodynamics import get_thermodynamic_information


from sens import NestedSamplingSAExact, HSASamplerCluster
from sens import get_all_normalmodes
from sens.models._lj_tools import LJClusterSENS

import _test_ns_lj

import _utils

class TestSENSExact_LJ(_test_ns_lj.TestNS_LJ):
    def setUp(self):
        self.seed = np.random.randint(1000000)
#        self.seed = 549670 # failed self.assertGreater(self.ns.count_sampled_minima, 0)
        print "seed", self.seed
        np.random.seed(self.seed)
        self.setUp1()
    
    def set_up_system(self):
        self.natoms = 6
        self.gmin = -12.7121
        self.system = LJClusterSENS(self.natoms, 2.5)
        self.ndof = 3*self.natoms - 6


    def setUp1(self, nproc=1):
        self.set_up_system()
        self.nreplicas = 10# * nproc
        self.stepsize = 0.01
        self.nproc = nproc
        
        self.database = self.system.create_database()
        # add some minima to the database
        bh = self.system.get_basinhopping(self.database, outstream=None)
        while self.database.number_of_minima() < 2:
            bh.run(1)
        # compute the thermodynamic information
        get_thermodynamic_information(self.system, self.database)
        get_all_normalmodes(self.system, self.database)
        

        self.minima = list(self.database.minima())
        assert self.database.number_of_minima() > 1, "%d minima" %  self.database.number_of_minima()
        
        self.mc_runner = self.system.get_mc_walker(mciter=200)

        self.energy_accuracy = 1e-4
        self.hsa_sampler = HSASamplerCluster(self.minima, self.system.k, copy_minima=True, 
                                             center_minima=True, 
                                             energy_accuracy=self.energy_accuracy, 
                                             compare_structures=self.system.get_compare_exact(), 
                                             mindist=self.system.get_mindist(), 
                                             minimizer=self.system.get_minimizer(), 
                                             debug=True)

        replicas = _utils.create_replicas(self.system, self.nreplicas)
        potential = self.system.get_potential()
        potential.get_energy = potential.getEnergy
        self.ns = NestedSamplingSAExact(replicas, self.mc_runner,
                                   self.hsa_sampler, potential,
                                   config_tests=self.system.get_config_tests(),
                                   nproc=nproc, verbose=True, iprint=1, debug=True)
        
        self.Emax0 = self.ns.replicas[-1].energy
        
        self.run_ns(max_iter=1000, Etol=.001)
    

    
    def test1(self):
        super(TestSENSExact_LJ, self).test1()
        self.assertGreater(self.ns.number_swaps_accepted(), 0)
        
        
        
        

class TestSENSExact_LJ_Par(TestSENSExact_LJ):
    def setUp(self):
        print "\n\nnproc ", 2
        self.setUp1(nproc=2)
    
    
if __name__ == "__main__":
    unittest.main()  
