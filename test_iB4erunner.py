#!/usr/bin/env python
import unittest
import os, shutil
import iB4erunner

class iB4erunnerTest(unittest.TestCase):
    def setUp(self):
        self.sagetarget = "test_result.sage"
        self.structdir = "test_results"
        self.paramdir = "Turner99"
        try:
            os.remove(self.sagetarget)
        except OSError:
            if os.path.isfile(self.sagetarget):
                pass

        try:
            os.makedirs(self.structdir)
        except OSError:
            if not os.path.isdir(self.structdir):
                raise

    def tearDown(self):
        try:
            shutil.rmtree(self.structdir)
        except OSError:
            if not os.path.isdir(self.structdir):
                raise
        try:
            os.remove(self.sagetarget)
        except OSError:
            if os.path.isfile(self.sagetarget):
                pass
    
    def test_combinatorial_sequence(self):
        seqfile = "test_data/test_tRNA.fasta"
        points = iB4erunner.run_iB4e(seqfile, self.sagetarget, self.paramdir, self.structdir)
        correct = [[8, 6, 24, 0], [8, 6, 24, 0], [0, 0, 0, 0], [0, 0, 0, 86.4], [2, 6, 6, -27.7], [8, 9, 25, 0], [0, 0, 0, 0], [1, 63, 3, 0], [8, 0, 24, 0], [8, 0, 24, 0], [1, 63, 3, 0], [2, 8, 6, -27.79999906794234]]
        self.assertEqual(correct, points)

    def test_arboreum_sequence(self):
        seqfile = "test_data/arboreum5S.fasta"
        points = iB4erunner.run_iB4e(seqfile, self.sagetarget, self.paramdir, self.structdir)
        correct = [[14, 10, 42, 0], [14, 10, 42, 0], [0, 0, 0, 0], [0, 0, 0, 141.38000363694223], [4, 20, 12, -45.14323869706544], [14, 17, 42, 0], [0, 0, 0, 0], [1, 108, 3, 0], [13, 0, 39, 0], [13, 0, 39, 0], [1, 108, 3, 0], [2, 21, 6, -49.70000088615559]]
        self.assertEqual(correct, points)

# Voodoo to make Python run the tests
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(iB4erunnerTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

