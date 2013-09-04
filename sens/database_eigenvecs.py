"""
attach eigenvalues and eigenvectors to the database
"""
import numpy as np

from sqlalchemy.orm import relationship, backref, deferred
from sqlalchemy import Column, Integer, Float, PickleType, ForeignKey

from pele.storage.database import Base
from pele.storage import Database


class HessianEigs(Base):
    """hessian eigenvalues and associated eigenvectors
        
    Parameters
    ----------
    m : Minimum object
    eigenvalues : array
        the list of eigenvalues
    eigenvectors : 2-d numpy array
        vectors[:,i] is the eigenvector associated to the normal mode frequency freq[i].
        the eigenvectors should be orthogonal and normalized
    """
    __tablename__ = "tbl_hessian_eigs"

    eigenvalues = deferred(Column(PickleType))
    eigenvectors = deferred(Column(PickleType))

    _minimum_id = Column(Integer, ForeignKey('tbl_minima._id'), primary_key=True)
    minimum = relationship("Minimum",
                            primaryjoin="Minimum._id==HessianEigs._minimum_id",
                            backref='hessian_eigs', uselist=False)
    def __init__(self, m, eigenvalues, eigenvectors):
        self.minimum = m
        self.eigenvalues = eigenvalues
        self.eigenvectors = np.copy(eigenvectors)


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
        
