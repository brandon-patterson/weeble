from sequence import AlignedSequence
from sequence_edit import SequenceReplacementEdit
from restriction_enzymes import RestrictionEnzyme
from wobble_cut_detector import WobbleCutDetector
import codons
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
            SequenceReplacementEdit(
                aligned, enzyme.sequence.reverse_complement(), 3)]
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

    def test_detect_cuts_handles_degenerate_bases(self):
        aligned = AlignedSequence('AGT')
        enzyme = RestrictionEnzyme('enzyme_x', 'NNN')
        actual_cuts = WobbleCutDetector().detect_cuts(aligned, enzyme)
        expected_cut_count = len([
            v for v in codons._encodings.values() if v == codons._encodings['AGT']
        ])
        self.assertEqual(len(actual_cuts), expected_cut_count)


if __name__ == '__main__':
    unittest.main()
