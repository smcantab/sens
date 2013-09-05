import unittest
import numpy as np

import sens.src.weighted_pick as wp


def weighted_pick(weights):
    """sample uniformly from the objects with given weights
    
    returns the selected index
    
    Parameters
    ----------
    
    weights : array of floats
        weight associated to each minimum, proportional to the volume of the basin
    """
    if len(weights) == 0:
        raise ValueError("weights must not have zero length")
    r = np.random.uniform(0., sum(weights))
    s = 0.0
#    print r, len(weights)    
    for i in range(len(weights)):
        s += weights[i]
        if r < s: return i
    return len(weights) - 1


class TestWeightedPick(unittest.TestCase):
    def setUp(self):
        self.weights = np.array(range(100), np.float)
        seed = np.random.randint(1)
        np.random.seed(seed)
        self.i1 = weighted_pick(self.weights)
        np.random.seed(seed)
        self.i2 = wp.weighted_pick_cython(self.weights)
    
    def test(self):
        self.assertEqual(self.i1, self.i2)
    
    
        
        
if __name__ == "__main__":
    unittest.main()  