"""
attach eigenvalues and eigenvectors to the database
"""
import numpy as np

from sqlalchemy.orm import relationship, backref, deferred
from sqlalchemy import Column, Integer, Float, PickleType, ForeignKey

from pele.storage.database import Base, Database, backref
from pele.thermodynamics import normalmodes as compute_normalmodes

class NormalModes(Base):
    """normal mode frequencies and associated eigenvectors
    
    actually frequencies are stored rather than eigenvalues.  Frequencies are
    eigenvalues with the appropriate weighting for (from the mass or from rigid bodies)
    
    Parameters
    ----------
    m : Minimum object
    freq : array
        the list of normal mode frequencies
    vectors : 2-d numpy array
        vectors[:,i] is the eigenvector associated to the normal mode frequency freq[i].
        the eigenvectors should be orthogonal and normalized
    nzero : integer, optional
        the number of zero frequencies (from global symmetries).  Default 0
    nnegative : integer, optional
        the number of negative frequencies (e.g. at a saddle point).  Default 0
    """
    __tablename__ = "tbl_normal_modes"

    freqs = deferred(Column(PickleType))
    vectors = deferred(Column(PickleType))
#    nzero = Column(Integer)
#    nnegative = Column(Integer)

    _minimum_id = Column(Integer, ForeignKey('tbl_minima._id'), primary_key=True)
    minimum = relationship("Minimum",
#                            primaryjoin="Minimum._id==NormalModes._minimum_id",
                            backref=backref('normal_modes', uselist=False))
    def __init__(self, m, freqs, vectors):
        self.minimum = m
        self.freqs = np.array(freqs)
        self.vectors = np.array(vectors)
        n = self.freqs.size
        assert self.vectors.shape == (n, n)
#        self.nzero = nzero
#        self.nnegative = nnegative


#class HessianEigs(Base):
#    """hessian eigenvalues and associated eigenvectors
#        
#    Parameters
#    ----------
#    m : Minimum object
#    eigenvalues : array
#        the list of eigenvalues
#    eigenvectors : 2-d numpy array
#        vectors[:,i] is the eigenvector associated to the normal mode frequency freq[i].
#        the eigenvectors should be orthogonal and normalized
#    """
#    __tablename__ = "tbl_hessian_eigs"
#
#    eigenvalues = deferred(Column(PickleType))
#    eigenvectors = deferred(Column(PickleType))
#
#    _minimum_id = Column(Integer, ForeignKey('tbl_minima._id'), primary_key=True)
#    minimum = relationship("Minimum",
#                            primaryjoin="Minimum._id==HessianEigs._minimum_id",
#                            backref='hessian_eigs', uselist=False)
#    def __init__(self, m, eigenvalues, eigenvectors):
#        self.minimum = m
#        self.eigenvalues = eigenvalues
#        self.eigenvectors = np.copy(eigenvectors)

def get_all_normalmodes(system, db):
    """
    for each minima in database, get all information necessary to compute the density of states
    
    Parameters
    ----------
    
    system : pele System
        is the particular system of interest, say LJCluster
    db : database of minima
    
    """

    # get the frequencies
#    print "getting the normal mode frequencies"
    for m in db.minima():
        # calculate the Hessian
        pot = system.get_potential()
        e, g, hess = pot.getEnergyGradientHessian(m.coords)
        
        # calculate the normal modes from the hessian
        freq, vectors = compute_normalmodes(hess, metric=None)
        nm = NormalModes(m, freq, vectors)
        db.session.add(nm)
        db.session.commit()

#        # calculate the eigenvalues and eigenvectors of the hessian and attach them to the database
#        eval, evec = get_eig(hess)
#        if len(m.hessian_eigs) == 0: 
#            nm = HessianEigs(m, eval, evec)

    db.session.commit()


if __name__ == "__main__":
    from pele.systems import LJCluster
    from pele.utils.hessian import get_sorted_eig
    system = LJCluster(20)
    db = system.create_database()#"test.sqlite")
    coords, energy = system.get_random_minimized_configuration()[:2]
    print energy
    m = db.addMinimum(energy, coords)
    
    pot = system.get_potential()
    e, g, hess = pot.getEnergyGradientHessian(coords)
    eval, evec = get_sorted_eig(hess)
    epair = NormalModes(m, eval, evec, nzero=6, nnegative=0)
    db.session.commit()
    print m.normal_modes[0].vectors.shape
    
    
    if False:
        print m.normal_modes 
        print m.normal_modes[0].freq, eval[0]
        epair = NormalModes(m, eval[1], evec[:,1])
        print m.normal_modes 
        print m.normal_modes[0].freq, eval[0]
        print m.normal_modes[1].freq, eval[1]
    elif False:
        for i in range(len(eval)):
            epair = NormalModes(m, eval[i], evec[:,i])
#            print m.normal_modes
        db.session.commit()
        print m.normal_modes[-1].freq
        print len(m.normal_modes)
        db.session.delete(m.normal_modes[0])
        db.session.commit()
        print len(m.normal_modes)

   