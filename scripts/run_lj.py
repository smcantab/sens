import argparse

from nested_sampling._nested_sampling import NestedSampling
from nested_sampling.nested_sampling_runner import run_nested_sampling

from sens.models._lj_tools import LJClusterSENS
from sens import NestedSamplingSA, NestedSamplingSAExact

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
    args = parser.parse_args()
    
    
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
    
    if args.sens_exact:
        
        ns = NestedSamplingSAExact(system, args.nreplicas, mcrunner,
                                   minima, energy_accuracy,
                                   **nskwargs)
    else:
        ns = NestedSampling(system, args.nreplicas, mcrunner, 
                            **nskwargs)
    
    run_nested_sampling(ns, label="lj"+str(args.natoms), etol=1e-3)
    

if __name__ == "__main__":
    main()