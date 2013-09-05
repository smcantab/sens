

#def get_normal_modes(system, db):
#    """
#    for each minima in database, get all information necessary to compute the density of states
#    
#    Parameters
#    ----------
#    
#    system : pele System
#        is the particular system of interest, say LJCluster
#    db : database of minima
#    
#    """
#    # get the point group information
#    print "getting the point group information"
#    determine_pgorder = system.get_pgorder()
#    for m in db.minima():
#        if m.pgorder is None:
#            m.pgorder = determine_pgorder(m.coords)
##        print m.pgorder
#    db.session.commit()
#
#    # get the frequencies
#    print "getting the normal mode frequencies"
#    for m in db.minima():
#         # fvibs calculates the eigenvalues and eigenvectors of the hessian and attach them to the database
#        if m.fvib is not None:
#            # assume we've already done this minimum 
#            #if len(m.hessian_eigs) > 0 and m.fvib is not None: 
#            continue 
#        # calculate the Hessian
#        pot = system.get_potential()
#        e, g, hess = pot.getEnergyGradientHessian(m.coords)
#        
#        # calculate the normal modes from the hessian
#        freq, evec = normalmodes(hess, metric=system.get_metric_tensor(m.coords))
#        # calculate the log product of the positive normal mode frequencies
#        n, lnf = logproduct_freq2(freq, system.nzero_modes)
#        m.fvib = lnf
#        
#        # calculate the eigenvalues and eigenvectors of the hessian and attach them to the database
#        eval, evec = get_eig(hess)
#        if len(m.hessian_eigs) == 0: 
#            nm = HessianEigs(m, eval, evec)
#
#    db.session.commit()
#    return db
