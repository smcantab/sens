import argparse

from nested_sampling import NestedSampling
from nested_sampling import run_nested_sampling, Replica

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
    parser.add_argument("--sens-approximate", action="store_true", help="use the approximate version of superposition enhanced nested sampling")
    parser.add_argument("--minprob", type=float, default=1., help="only for sens-approximate: probability of sampling from the database (by default <minprob>/K),"
                        "default value 1")
    parser.add_argument("--energy-offset", type=float, default=2.5, help="only for sens-approximate: value of Emax where the probabilty"
                        "starts to become non zero is <energy_max_database> + <energy_offset>")
    parser.add_argument("--db", type=str, help="location of the database", default="")
    parser.add_argument("--nminima", type=int, default=-1, help="number of minima from the database to use.  If negative, use all minima")
    parser.add_argument("--stop-crit", type=float, default=1e-3, help="run will terminate when stop_crit is larger than the difference between the maximum and minimum replica energies")
    parser.add_argument("--job-name", type=str, help="name of the job")
    parser.add_argument("--nsIP", type=str, help="IP address of the machine hosting the Name Server")
    args = parser.parse_args()
    print args
    
    if (args.sens_exact or args.sens_approximate) and args.db == "":
        raise Exception("for sens you must specify a database file")
    

    system = LJClusterSENS(args.natoms, args.radius)
    energy_accuracy = 1e-3

    if args.sens_exact or args.sens_approximate:
        database = system.create_database(args.db, createdb=False)
        if args.nminima <= 0 or args.nminima > database.number_of_minima(): 
            minima = database.minima()
        else:
            minima = database.minima()[:args.nminima]
        print "using", len(minima), "minima"
    
    
    mcrunner = system.get_mc_walker(args.mciter)
    nskwargs = dict(nproc=args.nproc, 
                    verbose=not args.q, iprint=args.iprint, job_name = args.job_name, nsIP = args.nsIP)
    
    # create the potential object
    potential = system.get_potential()
    potential.get_energy = potential.getEnergy

    # create the replicas
    replicas = []
    for i in xrange(args.nreplicas):
        x = system.get_random_configuration()
        e = potential.getEnergy(x)
        replicas.append(Replica(x, e))
        
        
    
    print "mciter", args.mciter
    print "radius", args.radius
    if args.sens_exact:
        hsa_sampler = HSASamplerCluster(minima, system.k, copy_minima=True, 
                center_minima=True, energy_accuracy=energy_accuracy, 
                compare_structures=system.get_compare_exact(), 
                mindist=system.get_mindist(), 
                minimizer=system.get_minimizer(), 
                debug=True)
        ns = NestedSamplingSAExact(replicas, mcrunner, hsa_sampler, potential,
                                   config_tests=system.get_config_tests(),
                                   **nskwargs)
    elif args.sens_approximate:
        ns = NestedSamplingSA(replicas, mcrunner, minima, system.k, 
                              config_tests=system.get_config_tests(),
                              minprob=args.minprob, energy_offset=args.energy_offset, 
                              copy_minima=True, center_minima=True, 
                              **nskwargs)
    else:
        ns = NestedSampling(replicas, mcrunner, 
                            **nskwargs)
    
    run_nested_sampling(ns, label="lj"+str(args.natoms), etol=args.stop_crit)
    

if __name__ == "__main__":
    main()
