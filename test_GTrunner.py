#!/usr/bin/env python
import unittest
import os
import GTrunner

class gtmfeTest(unittest.TestCase):
    def setUp(self):
        self.structtarget = "test_result.ct"
        self.paramdir = "Turner99"
        try:
            os.remove(self.structtarget)
        except OSError:
            if os.path.isfile(self.structtarget):
                pass

    def tearDown(self):
        try:
            os.remove(self.structtarget)
        except OSError:
            if os.path.isfile(self.structtarget):
                pass
    
    def test_combinatorial_sequence_classical(self):
        seqfile = "test_data/test_tRNA.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir)
        self.assertEqual([1, 10, 3], [result.multiloops, result.unpaired, result.branches])

    def test_combinatorial_sequence_variant(self):
        seqfile = "test_data/test_tRNA.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir, [-3.4, 0, -0.4, 1])
        self.assertEqual([2, 8, 6], [result.multiloops, result.unpaired, result.branches])

    def test_arboreum_classical(self):
        seqfile = "test_data/arboreum5S.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir)
        self.assertEqual([1, 25, 3], [result.multiloops, result.unpaired, result.branches])

    def test_arboreum_variant(self):
        seqfile = "test_data/arboreum5S.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir, [-3.4, 0, -0.4, 1])
        self.assertEqual([5, 23, 15], [result.multiloops, result.unpaired, result.branches])

# Voodoo to make Python run the tests
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(gtmfeTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
