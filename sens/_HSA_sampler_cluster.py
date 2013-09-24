import bisect
from collections import namedtuple
import numpy as np

from sens import SASampler
from twisted.python.formmethod import InputError


class _UnboundMinimumSmall(object):
    """represent a minimum object unbound from the database"""
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords


class _MinimaSearcher(object):
    """class to manage searching efficiently for minima in a database
    
    Parameters
    ----------
    minima : list of Mimumum objects
    energy_accuracy : float
        energy tolerance for when to consider that two minima are the same
    compare_structures : callable, `compare_structures(coords1, coords2)`
        a function that returns true if two structures are the same and false otherwise.
        This is used only if the two minima have energies that are within 
        energy_accuracy of each other.
    """
    def __init__(self, minima, energy_accuracy, compare_structures=None):
        self.minima = minima
        self.energy_accuracy = energy_accuracy
        self.compare_structures = compare_structures
        
        self.minima_dict = dict([(m.energy, m) for m in self.minima])
        self.energies = [m.energy for m in self.minima]
        self.energies.sort()
        
    def get_minimum(self, energy, coords):
        """return the minimum that matches energy and coords.  return None if there is no match"""
        # first get a list of minima energies that might be the same minima
        Emax = energy + self.energy_accuracy
        Emin = energy - self.energy_accuracy
        imax = bisect.bisect_right(self.energies, Emax)
        imin = bisect.bisect_left(self.energies, Emin)
        
#        print "found", imax - imin, "cadidates"
        RetVal = namedtuple("RetVal", ["minimum", "transformation"])
        retval = RetVal(minimum=None, transformation=None)
        mtest = _UnboundMinimumSmall(energy, coords)
        for e in self.energies[imin:imax]:
            assert Emin <= e <= Emax
            m = self.minima_dict[e]
            retval = RetVal(minimum=m, transformation=None)

            # compare the coordinates more carefully
            if self.compare_structures is None:
                break
            else:
                transform = self.compare_structures.find_transformation(m.coords, mtest.coords)
                if transform is not None:
                    # the structures match exactly.  
                    # transform is the transformation that turns coords into m.coords 
                    retval._replace(transformation=transform)
                    break
        
        return retval

class HSASamplerCluster(SASampler):
    """This class will manage sampling from the HSA as well as computing the HSA energy of a structure
    
    Notes
    -----
    Sampling from the HSA is done by the parent class SASampler.  This class is responsible for
    computing the energy of a given configuration within the HSA.  This can be seriously non-trivial
    because the full Hamiltonian often has symmetries that the HSA doesn't respect.  This class will
    handle one of the most difficult cases - a cluster with translational, rotational, permutational,
    and inversion degrees of freedom.  
    """
    def __init__(self, minima, k, copy_minima=True, center_minima=False,
                 energy_accuracy=1e-3, compare_structures=None, mindist=None,
                 minimizer=None, debug=True):
        super(HSASamplerCluster, self).__init__(minima, k, copy_minima=copy_minima, 
                                                center_minima=center_minima)
        
        if minimizer is None:
            raise ValueError("minimizer cannot be None")
        self.minima_searcher = _MinimaSearcher(minima, energy_accuracy, compare_structures=compare_structures)
        self.compare_structures = compare_structures
        self.mindist = mindist
        self.minimizer = minimizer
        self.debug = debug
        
        if self.debug:
            self._check_minima()

    def _check_minima(self):
        for m in self.minima:
            assert m.energy is not None
            assert m.fvib is not None
            assert m.pgorder is not None
            assert m.normal_modes is not None

 
    def compute_energy_in_HSA(self, energy, x):
        """compute the energy of the coordinates in replica.coords in the harmonic superposition approximation
        
        This is where most of the difficulty is in the algorithm.  
        
        This is also where all the system dependence is
        
        Returns
        -------
        energy : the energy of x in the HSA.  None if no energy can be found (if the closest minimum is not in the database).
        """
        x = x.copy()
        # quench to nearest minimum
        qresult = self.minimizer(x)
        
        # check if that minimum is in the database.  reject if not
        m, transformation = self.minima_searcher.get_minimum(qresult.energy, qresult.coords)
        if m is None:
            print "rejecting.  minima not in database"
            return None
        
        # put replica.coords into best alignment with the structure stored in m.coords
        # this involves accounting for symmetries of the Hamiltonian like translational, 
        # rotational and permutations symmetries.  You can use the coordinates in qresult.coords
        # to help find the best permutation.  Ultimately it must be aligned with the structure in m.coords
        # The hessian eigenvectors were computed with a given permutation
        # e.g. account for trivial translational and rotational degrees of freedom
        if transformation is not None:
            # transformation is the set of transformations that put qresult.coords into exact
            # alignment with m.coords.  If we apply these transformations to replica.x then
            # replica.x will be in good (although not perfect) alignment already with m.coords 
            x = self.compare_structures.apply_transformation(x, transformation)
            
        # Do a final round of optimization to further improve the alignment
        if self.mindist is not None:
            xcopy = x.copy()
            dist, x0, x = self.mindist(m.coords.copy(), x)
            if self.debug:
                diff = np.linalg.norm(x0 - m.coords)
                if diff > .01:
                    with open("error.xyz", "w") as fout:
                        from pele.utils.xyz import write_xyz
                        write_xyz(fout, m.coords, title="x1 initial")
                        write_xyz(fout, x0, title="x1 final")
                        write_xyz(fout, xcopy, title="x2 initial")
                        write_xyz(fout, x, title="x2 final")
                        
                    raise Exception("warning, mindist appears to have changed x0.  the norm of the difference is %g" % diff)
                    
                assert (x0 == m.coords).all()
        
        energy = self.compute_energy(x, m)
        
        return energy    
