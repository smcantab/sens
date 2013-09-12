from sens import get_all_normalmodes
from pele.thermodynamics import get_thermodynamic_information

def build_database(system, nminima, dbfname=None):
    if dbfname is None:
        db = system.create_database()
    else:
        db = system.create_database(dbfname)
    
    bh = system.get_basinhopping(db)
    
    while db.number_of_minima() < nminima:
        bh.run(1)
    
    get_thermodynamic_information(system, db)
    get_all_normalmodes(system, db)
    
    return db
        