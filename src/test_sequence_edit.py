from sequence import AlignedSequence
from sequence import Sequence
from sequence_edit import SequenceReplacementEdit
import unittest


class TestSequenceEdit(unittest.TestCase):
    def test_init_enforces_lengths(self):
        aligned_sequence = AlignedSequence('AAACCCGGGTTT')
        override_sequence = Sequence('AAAA')
        self.assertRaises(AssertionError, lambda: SequenceReplacementEdit(
            aligned_sequence, override_sequence, 9))

    def test_eq(self):
        edit = SequenceReplacementEdit(AlignedSequence('ACT'), Sequence('A'), 2)
        same = SequenceReplacementEdit(AlignedSequence('ACT'), Sequence('A'), 2)
        diff = SequenceReplacementEdit(AlignedSequence('ACT'), Sequence('G'), 2)
        self.assertEqual(edit, same)
        self.assertNotEqual(edit, diff)

    def test_str(self):
        aligned_sequence = AlignedSequence('AAACCCGGGTTT')
        override_sequence = Sequence('AA')
        edit = SequenceReplacementEdit(aligned_sequence, override_sequence, 5)
        self.assertEqual(str(edit),
                         'edit index: 4\tCCCGGG -> CCaaGG\tbases changed: 2\tabs usage shift: inf')

    def test_get_number_of_bases_modified(self):
        aligned_sequence = AlignedSequence('AAACCCGGGTTT')
        override_sequence = Sequence('AA')
        edit = SequenceReplacementEdit(aligned_sequence, override_sequence, 5)
        self.assertEqual(edit.get_number_of_bases_modified(), 2)

    def test_get_number_of_amino_modified(self):
        aligned_sequence = AlignedSequence('AAACCCGGGTTT')
        override_sequence = Sequence('AA')
        edit = SequenceReplacementEdit(aligned_sequence, override_sequence, 5)
        # CCC and CCa both code P
        self.assertEqual(edit.get_number_of_aminos_modified(), 1)

    def test_get_abs_usage_shift(self):
        aat = AlignedSequence('aat')
        aac = AlignedSequence('aac')
        aataac = AlignedSequence('aataac')
        aacaat = AlignedSequence('aacaat')
        self.assertEqual(SequenceReplacementEdit(aat, aat, 0)
                         .get_abs_usage_shift(), 0)
        self.assertAlmostEqual(SequenceReplacementEdit(aat, aac, 0)
                               .get_abs_usage_shift(), 0.29)
        self.assertAlmostEqual(SequenceReplacementEdit(aac, aat, 0)
                               .get_abs_usage_shift(), 0.29)
        self.assertAlmostEqual(SequenceReplacementEdit(aacaat, aataac, 0)
                               .get_abs_usage_shift(), 0.58)


if __name__ == '__main__':
    unittest.main()
