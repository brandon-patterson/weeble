"""
A utility that searches for possible restriction enzyme cut-sites in a DNA
sequence, allowing for base substitutions (provided that the resulting amino
acid chain remains unaltered.
"""

from sequence import Sequence
from wobble_cut_detector import WobbleCutDetector
import codons
import os
import sys
import restriction_enzymes as enzymes

_default_args = {
    'sequence': 'sequence',
    'usage': 'human',
    'encodings': 'encodings',
    'restriction_enzymes': 'restriction_enzymes'
}


def _set_pwd_to_main():
    """Sets pwd to this file's location to ensure relative config paths work."""
    main_path = os.path.realpath(sys.argv[0])
    os.chdir(os.path.dirname(main_path))


def _get_inputs():
    arg_dict = _default_args
    args = sys.argv[1:]

    for i in range(len(args)):
        arg = args[i]

        # allow the first arg to be an unlabeled literal sequence
        if i == 0 and '=' not in arg:
            arg_dict['literal_sequence'] = args[0]
        else:
            k, v = arg.split('=')
            if k not in _default_args.keys():
                raise KeyError('argument "{}" not supported'.format(k))
            arg_dict[k] = v
    return arg_dict


def _load_sequence(sequence_file):
    seq = ''
    with open('../configs/{}.txt'.format(sequence_file), 'r') as infile:
        for line in infile:
            seq += line
    return Sequence(seq).align()


if __name__ == '__main__':
    _set_pwd_to_main()
    inputs = _get_inputs()

    codons.load_encodings(inputs['encodings'])
    codons.load_usage(inputs['usage'])
    aligned_seq = None
    if 'literal_sequence' in inputs.keys():
        aligned_seq = Sequence(inputs['literal_sequence']).align()
    else:
        aligned_seq = _load_sequence(inputs['sequence'])

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
