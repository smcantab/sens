import argparse

import sens._database_normalmodes
from sens._database_normalmodes import NormalModes

from pele.systems import LJCluster
from sens import get_all_normalmodes
from pele.thermodynamics import get_thermodynamic_information



def main():
    parser = argparse.ArgumentParser(description="do nested sampling on a Lennard Jones cluster")
    parser.add_argument("natoms", type=int, help="number of atoms")
    parser.add_argument("db", type=str, help="location of the database")
    parser.add_argument("--normalmodes", type=str, help="store the full normalmodes as well")
    parser.add_argument("--recalculate", action="store_true", help="recompute the thermodynamic information")
    args = parser.parse_args()
    
    

    system = LJCluster(args.natoms,)

    if False:
        # build the database    
        database = system.create_database(args.db, createdb=True)
        bh = system.get_basinhopping(database)
        bh.run(1000)
        
    
    else:
        database = system.create_database(args.db, createdb=False)

    assert database.number_of_minima() > 0
    
    print "computing the vibrational free energy and the point group order"
    get_thermodynamic_information(system, database, recalculate=args.recalculate)
    print "computing the normal modes"
    get_all_normalmodes(system, database)

    
    

if __name__ == "__main__":
    main()