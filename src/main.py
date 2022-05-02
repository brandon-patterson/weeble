"""
A utility that searches for possible restriction enzyme cut-sites in a DNA
sequence, allowing for base substitutions (provided that the resulting amino
acid chain remains unaltered.
"""

from sequence import Sequence
from wobble_cut_detector import WobbleCutDetector
import codons
import restriction_enzymes as enzymes

# Must be properly left-aligned. (Use leading underscores to align first codon.)
sequence_of_interest = \
    'AGTGGCTCCTAGAGGCTCGAATCAGCTATCAGCCGTACGGCATCATCAAACTTCTGGGGCTGC'

if __name__ == '__main__':
    aligned_seq = Sequence(sequence_of_interest).align()
    print('Original base sequence:\n{}'.format(aligned_seq))
    print()
    print('Original amino chain: \n\t{}'.format(aligned_seq.get_amino_string()))
    print()
    print('Checking sequence for potential cuts (leaving aminos unchanged)...')
    print('(usage table: {})'.format(codons._usage_source))
    print()

    wobbler = WobbleCutDetector()

    for enzyme in enzymes.get_all_enzymes():
        wobble_cuts = wobbler.detect_cuts(aligned_seq, enzyme)
        wobble_cuts.sort(key=lambda cut: (
            cut.get_number_of_bases_modified(), cut.get_abs_usage_shift()))
        if len(wobble_cuts) > 0:
            print('possible cuts for {}:'.format(enzyme.name))
            for cut_edit in wobble_cuts:
                print('  {}'.format(cut_edit))

            print()
