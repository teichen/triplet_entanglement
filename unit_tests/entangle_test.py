import numpy as np
import unittest

class EntanglementTester(unittest.TestCase):
    """ test entanglement calculation """

    def test_entanglement(self):
        """
        """
        pass

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(EntanglementTester)
    unittest.TextTestRunner().run(suite)
