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
        correct = [[8, 6, 24, 0], [8, 6, 24, 0], [0, 0, 0, 0], [2, 1, 11, 86.5], [1, 23, 3, -25.1], [1, 63, 3, 5.083333333333333], [8, 0, 24, 52.0], [1, 3, 13, 59.4], [0, 0, 0, 86.4], [0, 0, 0, 86.4]]
        self.assertEqual(correct, points)

    def test_arboreum_sequence(self):
        seqfile = "test_data/arboreum5S.fasta"
        points = iB4erunner.run_iB4e(seqfile, self.sagetarget, self.paramdir, self.structdir)
        correct = [[14, 10, 42, 0], [14, 10, 42, 0], [0, 0, 0, 0], [2, 0, 11, 141.38000363694223], [3, 38, 9, -47.5], [1, 108, 3, 14.4], [13, 0, 39, 69.04002403846154], [1, 12, 19, 71.3], [1, 0, 3, 141.37999334221038], [1, 0, 3, 141.37999334221038]]
        self.assertEqual(correct, points)

# Voodoo to make Python run the tests
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(iB4erunnerTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

