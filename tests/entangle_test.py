import sys
sys.path.insert(1, '../')
import numpy as np
import unittest
from MMax import MMax

class EntanglementTester(unittest.TestCase):
    """ test entanglement calculation """

    def test_entanglement(self):
        """
        """
        n  = 101 # size of the lattice
        nr = 20
        dr = 2.0 / (nr-1) # input triplet separation resolution
        r  = np.zeros(nr)

        for idx in range(nr):
            r[idx] = idx * dr

        ymean = MMax(r, n)

        self.assertAlmostEqual(ymean[0],  0.008163, places=4)
        self.assertAlmostEqual(ymean[2],  0.010204, places=4)
        self.assertAlmostEqual(ymean[-2], 0.226531, places=4)
        self.assertAlmostEqual(ymean[-1], 0.202041, places=4)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(EntanglementTester)
    unittest.TextTestRunner().run(suite)
