from sequence import Sequence
from wobble_cut_detector import WobbleCutDetector
import restriction_enzymes as enzymes


# Must be properly aligned. (Use leading underscores to align first codon.)
sequence_of_interest = 'AGTGGCTCCTAGAGGCTCGAATCAGCTATC'


if __name__ == '__main__':
    aligned_seq = Sequence(sequence_of_interest).align()
    print('original sequence: {}'.format(aligned_seq.bases))
    print('original amino chain: {}'.format(aligned_seq.get_amino_string()))
    print()
    print('Checking sequence for potential cuts (leaving aminos unchanged)...')
    print()

    wobbler = WobbleCutDetector()

    for enzyme in enzymes.get_all_enzymes():
        wobble_cuts = wobbler.detect_cuts(aligned_seq, enzyme)
        if len(wobble_cuts) > 0:
            print('possible cuts for {}:'.format(enzyme.name))
            for cut_edit in wobble_cuts:
                print('  {}'.format(cut_edit))

            print()
