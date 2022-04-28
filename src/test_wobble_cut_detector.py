from sequence import AlignedSequence
from sequence_edit import SequenceReplacementEdit
from restriction_enzymes import RestrictionEnzyme
from wobble_cut_detector import WobbleCutDetector
import unittest


class TestWobbleCutDetector(unittest.TestCase):
    def test_detect_cuts_no_match(self):
        aligned = AlignedSequence('AAATTT')
        enzyme = RestrictionEnzyme('enzyme_x', 'GGG')
        self.assertEqual(WobbleCutDetector().detect_cuts(aligned, enzyme), [])

    def test_detect_cuts_one_match(self):
        aligned = AlignedSequence('AGCGGCTTT')
        enzyme = RestrictionEnzyme('enzyme_x', 'GGG')
        actual_cuts = WobbleCutDetector().detect_cuts(aligned, enzyme)
        expected_cuts = [SequenceReplacementEdit(aligned, enzyme.sequence, 3)]
        self.assertEqual(actual_cuts, expected_cuts)

    def test_detect_cuts_one_reverse_match(self):
        aligned = AlignedSequence('AGCGGCTTT')
        enzyme = RestrictionEnzyme('enzyme_x', 'CCC')
        actual_cuts = WobbleCutDetector().detect_cuts(aligned, enzyme)
        expected_cuts = [
            SequenceReplacementEdit(aligned, enzyme.sequence.invert(), 3)]
        self.assertEqual(actual_cuts, expected_cuts)

    def test_detect_cuts_two_matches(self):
        aligned = AlignedSequence('AAAGGGTTT')
        enzyme = RestrictionEnzyme('enzyme_x', 'GGG')
        actual_cuts = WobbleCutDetector().detect_cuts(aligned, enzyme)
        expected_cuts = [
            SequenceReplacementEdit(aligned, enzyme.sequence, 2),
            SequenceReplacementEdit(aligned, enzyme.sequence, 3),
        ]
        self.assertEqual(actual_cuts, expected_cuts)

    def test_detect_cuts_palindromes_counted_once(self):
        aligned = AlignedSequence('AAATTT')
        enzyme = RestrictionEnzyme('enzyme_x', 'AT')
        actual_cuts = WobbleCutDetector().detect_cuts(aligned, enzyme)
        expected_cuts = [SequenceReplacementEdit(aligned, enzyme.sequence, 2)]
        self.assertEqual(actual_cuts, expected_cuts)


if __name__ == '__main__':
    unittest.main()
