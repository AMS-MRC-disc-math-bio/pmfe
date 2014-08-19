#!/usr/bin/env python
import unittest
import os
import GTrunner
from fractions import Fraction
from gtmfe import gtmfe

class gtmfeTest(unittest.TestCase):
    def setUp(self):
        # The default vector is [3.4, 0, 0.4, 1]. This is a weird one.
        self.weirdvector = gtmfe.ParameterVector()
        self.weirdvector.set_from_fractions(Fraction(18,3), Fraction(-12,5), Fraction(8,7), Fraction(-2,3))
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

    def parse_result(self, result):
        resultdict = result.get_python_numbers()
        multiloops = resultdict["multiloops"]
        unpaired = resultdict["unpaired"]
        branches = resultdict["branches"]
        return [multiloops, unpaired, branches]        
    
    def test_combinatorial_sequence_classical(self):
        seqfile = "test_data/test_tRNA.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir)
        self.assertEqual([1, 10, 3], self.parse_result(result))

    def test_combinatorial_sequence_variant(self):
        seqfile = "test_data/test_tRNA.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir, self.weirdvector)
        self.assertEqual([2, 8, 6], self.parse_result(result))

    def test_arboreum_classical(self):
        seqfile = "test_data/arboreum5S.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir)
        self.assertEqual([1, 25, 3], self.parse_result(result))

    def test_arboreum_variant(self):
        seqfile = "test_data/arboreum5S.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir, self.weirdvector)
        self.assertEqual([5, 23, 15], self.parse_result(result))

# Voodoo to make Python run the tests
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(gtmfeTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
