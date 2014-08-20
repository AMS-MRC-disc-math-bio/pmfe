#!/usr/bin/env python
import unittest
import os
import GTrunner
from fractions import Fraction
from gtmfe import gtmfe

class gtmfeTest(unittest.TestCase):
    def setUp(self):
        # The default vector is [3.4, 0, 0.4, 1]. Here are some others.
        # Vectors can be expressed with decimal strings
        self.variant_negative = gtmfe.ParameterVector()
        self.variant_negative.set_from_words('-3.4', '0', '-0.4', '1')

        # Vectors can also be expressed as exact rationals
        self.variant_madeup = gtmfe.ParameterVector()
        self.variant_madeup.set_from_fractions(Fraction(18,5), Fraction(-12,5), Fraction(8,7), Fraction(-2,3))
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

    def test_combinatorial_sequence_variant_negative(self):
        seqfile = "test_data/test_tRNA.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir, self.variant_negative)
        self.assertEqual([2, 8, 6], self.parse_result(result))

    def test_combinatorial_sequence_variant_madeup(self):
        seqfile = "test_data/test_tRNA.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir, self.variant_madeup)
        self.assertEqual([1, 63, 3], self.parse_result(result))

    def test_arboreum_classical(self):
        seqfile = "test_data/arboreum5S.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir)
        self.assertEqual([1, 25, 3], self.parse_result(result))

    def test_arboreum_variant_negative(self):
        seqfile = "test_data/arboreum5S.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir, self.variant_negative)
        self.assertEqual([5, 23, 15], self.parse_result(result))

    def test_arboreum_variant_madeup(self):
        seqfile = "test_data/arboreum5S.fasta"
        result = GTrunner.run_gtmfe(seqfile, self.structtarget, self.paramdir, self.variant_madeup)
        self.assertEqual([1, 108, 3], self.parse_result(result))

# Voodoo to make Python run the tests
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(gtmfeTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

