"""
this test will use the harmonic superposition approximation generate random configurations for the nested sampling algorithm.  
This should be exactly equal to the directly computable HSA results 
"""
from pele.optimize import Result

from nested_sampling import NestedSampling, run_nested_sampling, Replica

from _utils import build_database
from sens.models import LJClusterSENS
from sens import SASampler

    

class SamplerSystem(object):
    def __init__(self, natoms=6):
        self.system = LJClusterSENS(natoms, 5.)
        self.ndof = natoms * 3 - 6
        self.natoms = natoms
        
        self.database = build_database(self.system, nminima=2)
        
        self.sampler = SASampler(self.database.minima(), self.ndof)
    
    def sample_configuration(self, Emax):
        m, x = self.sampler.sample_coords(Emax)
        E = self.sampler.compute_energy(x, m)
        return E, x
    
    def __call__(self, *args, **kwargs):
        Emax = args[2]
#        print "walking with Emax", Emax
        E, x = self.sample_configuration(Emax)
        
        res = Result(x=x, energy=E, naccept=5, nsteps=10, Emax=Emax)
        return res
        
    
    

def main():
    system = SamplerSystem()
    
    
    nreplicas = 10000
    replicas = [Replica(x, energy) for energy, x in [system.sample_configuration(10.) for i in range(nreplicas)]]
    
    mc_walker = system 
    ns = NestedSampling(system, nreplicas=nreplicas, mc_walker=system, nproc=1, verbose=True, iprint=100, replicas=replicas)
    run_nested_sampling(ns, label="test")

if __name__ == "__main__":
    main()
        