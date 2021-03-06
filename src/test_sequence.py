from codons import Codon
from sequence import AlignedSequence
from sequence import Sequence
import unittest


class TestSequence(unittest.TestCase):
    def test_init_validates_bases(self):
        self.assertRaises(AssertionError, lambda: Sequence('XYZ'))
        self.assertRaises(AssertionError, lambda: Sequence('WIFI'))
        # Valid
        Sequence('_ACTUCGRYKMBVDHGSWN')

    def test_init_is_case_insensitive(self):
        self.assertEqual(Sequence('AGTC'), Sequence('agtc'))

    def test_eq(self):
        self.assertEqual(Sequence('ACT'), Sequence('ACT'))
        self.assertNotEqual(Sequence('ACT'), Sequence('CAT'))
        self.assertNotEqual(Sequence('AAA'), Sequence('AAAA'))

    def test_hash(self):
        self.assertEqual(hash(Sequence('ACT')), hash(Sequence('ACT')))
        self.assertNotEqual(hash(Sequence('ACT')), hash(Sequence('CAT')))
        sequence_set = set()
        sequence_set.add(Sequence('ACT'))
        sequence_set.add(Sequence('CAT'))
        sequence_set.add(Sequence('CAT'))
        sequence_set.add(Sequence('CAN'))
        self.assertEqual(len(sequence_set), 3)



    def test_len(self):
        self.assertEqual(len(Sequence('ACT')), 3)
        self.assertEqual(len(Sequence('CATG')), 4)

    def test_str(self):
        self.assertEqual(str(Sequence('ACT')), 'ACT')

    def test_reverse_complement(self):
        self.assertEqual(Sequence('ACTG_').reverse_complement(),
                         Sequence('_CAGT'))
        # Non-standard bases can also be inverted.
        self.assertEqual(Sequence('BDHKMNRSVWY').reverse_complement(),
                         Sequence('RWBSYNKMDHV'))

    def test_align_pads_missing_bases(self):
        self.assertEqual(Sequence('ACTA').align(), AlignedSequence('ACTA__'))

    def test_is_degenerate(self):
        self.assertFalse(Sequence('ACGTU_').is_degenerate())
        for c in 'BDHKMNRSVWY':
            self.assertTrue(Sequence(c).is_degenerate())

    def test_get_primitive_sequences(self):
        self.assertListEqual(Sequence('ACT').get_primitive_sequences(),
                             [Sequence('ACT')])
        self.assertListEqual(Sequence('ANT').get_primitive_sequences(),
                             [Sequence('AAT'), Sequence('ACT'),
                              Sequence('AGT'), Sequence('ATT')])
        self.assertListEqual(Sequence('RAM').get_primitive_sequences(),
                             [Sequence('AAA'), Sequence('AAC'),
                              Sequence('GAA'), Sequence('GAC')])


class TestAlignedSequence(unittest.TestCase):
    def test_init_validates_length(self):
        self.assertRaises(AssertionError, lambda: AlignedSequence('AA'))
        self.assertRaises(AssertionError, lambda: AlignedSequence('AAAA'))

    def test_init_populates_codons(self):
        seq = AlignedSequence('ACTGGC')
        expected_codons = [Codon('ACT'), Codon('GGC')]
        self.assertEqual(seq.codons, expected_codons)

    def test_str(self):
        long_seq = AlignedSequence('A' * 99)
        expected = \
            '1\tAAA AAA AAA AAA AAA : AAA AAA AAA AAA AAA\n' \
            + '31\tAAA AAA AAA AAA AAA : AAA AAA AAA AAA AAA\n' \
            + '61\tAAA AAA AAA AAA AAA : AAA AAA AAA AAA AAA\n' \
            + '91\tAAA AAA AAA'

        self.assertEqual(str(long_seq), expected)

    def test_get_amino_string(self):
        self.assertEqual(
            AlignedSequence('ACTGGC' * 6 + 'TT_').get_amino_string(),
            'TGTGT : GTGTG : TG?')

    def test_encodes_same_aminos(self):
        seq = AlignedSequence('ACTGGCTT_')
        matching = AlignedSequence('ACGGGATT_')
        not_matching = AlignedSequence('AAAGGGCCC')
        self.assertTrue(seq.encodes_same_aminos(matching))
        self.assertFalse(seq.encodes_same_aminos(not_matching))


if __name__ == '__main__':
    unittest.main()
