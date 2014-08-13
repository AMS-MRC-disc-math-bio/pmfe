#!/usr/bin/env python
import unittest
import os
import RNAscorer

class RNAscorerTest(unittest.TestCase):
    def test_combinatorial_structure(self):
        structfile = "test_data/test_tRNA.correct.ct"
        scores = RNAscorer.score_file(structfile)
        self.assertEqual((1, 10, 3), scores)

    def test_arboreum_structure(self):
        structfile = "test_data/arboreum5S.correct.ct"
        scores = RNAscorer.score_file(structfile)
        self.assertEqual((1, 25, 3), scores)

# Voodoo to make Python run the tests
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(RNAscorerTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
