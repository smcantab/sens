import argparse

from nested_sampling._nested_sampling import NestedSampling
from nested_sampling import run_nested_sampling

from sens.models._lj_tools import LJClusterSENS
from sens import NestedSamplingSA, NestedSamplingSAExact, HSASamplerCluster

def main():
    parser = argparse.ArgumentParser(description="do nested sampling on a Lennard Jones cluster")
    parser.add_argument("natoms", type=int, help="number of atoms")
    parser.add_argument("-K", "--nreplicas", type=int, help="number of replicas", default=300)
    parser.add_argument("-n", "--mciter", type=int, default=1000, help="number of steps in the monte carlo walk")
    parser.add_argument("-P", "--nproc", type=int, help="number of processors", default=1)
    parser.add_argument("-q", action="store_true", help="turn off verbose printing of information at every step")
    parser.add_argument("--radius", type=float, default=2.5, help="maintain atoms in a sphere of this radius")
    parser.add_argument("--iprint", type=int, default=1, help="if verbose, status messages will be printed every iprint steps")
    parser.add_argument("--sens-exact", action="store_true", help="use the exact version of superposition enhanced nested sampling")
    parser.add_argument("--db", type=str, help="location of the database", default="")
    parser.add_argument("--nminima", type=int, default=-1, help="number of minima from the database to use.  If negative, use all minima")
    parser.add_argument("--stop-crit", type=float, default=1e-5, help="run will terminate when stop_crit is larger than the difference between the maximum and minimum replica energies")
    args = parser.parse_args()
    print args
    
    if args.sens_exact and args.db == "":
        raise Exception("for sens you must specify a database file")
    

    system = LJClusterSENS(args.natoms, args.radius)
    energy_accuracy = 1e-3

    if args.sens_exact:
        database = system.create_database(args.db, createdb=False)
        if args.nminima <= 0 or args.nminima > database.number_of_minima(): 
            minima = database.minima()
        else:
            minima = database.minima()[:args.nminima]

    
    
    mcrunner = system.get_mc_walker(args.mciter)
    nskwargs = dict(nproc=args.nproc, 
                    verbose=not args.q, iprint=args.iprint)
    
    print "mciter", args.mciter
    print "radius", args.radius
    if args.sens_exact:
        hsa_sampler = HSASamplerCluster(minima, system.k, copy_minima=True, 
                center_minima=True, energy_accuracy=energy_accuracy, 
                compare_structures=system.get_compare_exact(), 
                mindist=system.get_mindist(), 
                minimizer=system.get_minimizer(), 
                debug=True)
        ns = NestedSamplingSAExact(system, args.nreplicas, mcrunner, hsa_sampler,
                                   config_tests=system.get_config_tests(),
                                   **nskwargs)
    else:
        ns = NestedSampling(system, args.nreplicas, mcrunner, 
                            **nskwargs)
    
    run_nested_sampling(ns, label="lj"+str(args.natoms), etol=args.stop_crit)
    

if __name__ == "__main__":
    main()