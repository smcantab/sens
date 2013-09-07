import unittest
import numpy as np

from pele.thermodynamics import get_thermodynamic_information

from sens._sens_exact import NestedSamplingSAExact
from sens import get_all_normalmodes
from sens.models._lj_tools import LJClusterSENS

import _test_ns_lj

class TestSENS_LJ(_test_ns_lj.TestNS_LJ):
    def setUp(self):
        self.setUp1()
    
    def set_up_system(self):
        self.natoms = 6
        self.gmin = -12.7121
        self.system = LJClusterSENS(self.natoms, 2.5)


    def setUp1(self, nproc=1):
        self.set_up_system()
        self.nreplicas = 10
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
        
        self.mc_runner = self.system.get_mc_walker(mciter=100)

        self.energy_accuracy = 1e-4
        self.ns = NestedSamplingSAExact(self.system, self.nreplicas, self.mc_runner,
                                   self.minima, self.energy_accuracy, 
                                   compare_minima=self.system.get_compare_minima(), 
                                   mindist=self.system.get_mindist(),
                                   config_tests = self.system.get_config_tests(),
                                   stepsize=0.1, nproc=nproc, verbose=True)
        
        self.Emax0 = self.ns.replicas[-1].energy
        
        self.run_ns(max_iter=1000, Etol=.01)
    
#    def test2(self):
#        self.assertGreater(self.ns.count_sampled_minima, 0)
    

#class testSENS_LJ_Par(TestSENS_LJ):
#    def setUp(self):
#        self.setUp1(nproc=3)
    
    
if __name__ == "__main__":
    unittest.main()  
