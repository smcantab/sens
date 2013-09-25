from sens import get_all_normalmodes
from pele.thermodynamics import get_thermodynamic_information
from nested_sampling import Replica

def build_database(system, nminima, dbfname=None, maxiter=1000):
    if dbfname is None:
        db = system.create_database()
    else:
        db = system.create_database(dbfname)
    
    bh = system.get_basinhopping(db, outstream=None)
    
    i = 0
    while db.number_of_minima() < nminima:
        bh.run(1)
        i += 1
        if i >= maxiter:
            break
    
    get_thermodynamic_information(system, db)
    get_all_normalmodes(system, db)
    
    return db

def create_replicas(system, nreplicas):
    replicas = []
    for i in xrange(nreplicas):
        x = system.get_random_configuration()
        replicas.append(Replica(x, system.get_energy(x)))
    return replicas
