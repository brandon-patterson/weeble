from codons import Codon
from sequence import AlignedSequence
from sequence import Sequence
import unittest

class TestSequence(unittest.TestCase):
    def test_init_validates_bases(self):
        self.assertRaises(AssertionError, lambda: Sequence('XYZ'))
        self.assertRaises(AssertionError, lambda: Sequence('BBB'))
        # Valid
        Sequence('_ACTG')

    def test_init_is_case_insensitive(self):
        self.assertEqual(Sequence('AGTC'), Sequence('agtc'))

    def test_eq(self):
        self.assertEqual(Sequence('ACT'), Sequence('ACT'))
        self.assertNotEqual(Sequence('ACT'), Sequence('CAT'))
        self.assertNotEqual(Sequence('AAA'), Sequence('AAAA'))

    def test_len(self):
        self.assertEqual(len(Sequence('ACT')), 3)
        self.assertEqual(len(Sequence('CATG')), 4)

    def test_invert(self):
        self.assertEqual(Sequence('ACTG_').invert(), Sequence('_CAGT'))

    def test_align_pads_missing_bases(self):
        self.assertEqual(Sequence('ACTA').align(), AlignedSequence('ACTA__'))


class TestAlignedSequence(unittest.TestCase):
    def test_init_validates_length(self):
        self.assertRaises(AssertionError, lambda: AlignedSequence('AA'))
        self.assertRaises(AssertionError, lambda: AlignedSequence('AAAA'))

    def test_init_populates_codons(self):
        seq = AlignedSequence('ACTGGC')
        expected_codons = [Codon('ACT'), Codon('GGC')]
        self.assertEqual(seq.codons, expected_codons)

    def test_get_amino_string(self):
        self.assertEqual(AlignedSequence('ACTGGCTT_').get_amino_string(), 'TG?')

    def test_encodes_same_aminos(self):
        seq = AlignedSequence('ACTGGCTT_')
        matching = AlignedSequence('ACGGGATT_')
        not_matching = AlignedSequence('AAAGGGCCC')
        self.assertTrue(seq.encodes_same_aminos(matching))
        self.assertFalse(seq.encodes_same_aminos(not_matching))


if __name__ == '__main__':
    unittest.main()
