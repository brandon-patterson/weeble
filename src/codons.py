import base_utils
import itertools

_UNKNOWN_AMINO = '?'

_encodings = {
    'AAA': 'K',
    'AAC': 'N',
    'AAG': 'K',
    'AAT': 'N',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'AGA': 'R',
    'AGC': 'S',
    'AGG': 'R',
    'AGT': 'S',
    'ATA': 'I',
    'ATC': 'I',
    'ATG': 'M',
    'ATT': 'I',
    'CAA': 'Q',
    'CAC': 'H',
    'CAG': 'Q',
    'CAT': 'H',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'GAA': 'E',
    'GAC': 'D',
    'GAG': 'E',
    'GAT': 'D',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'TAA': 'Ochre',  # Use a matching symbol to make stops interchangeable
    'TAC': 'Y',
    'TAG': 'Amber',  # Use a matching symbol to make stops interchangeable
    'TAT': 'Y',
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TGA': 'Opal',  # Use a matching symbol to make stops interchangeable
    'TGC': 'C',
    'TGG': 'W',
    'TGT': 'C',
    'TTA': 'L',
    'TTC': 'F',
    'TTG': 'L',
    'TTT': 'F',
}

_usage_source = 'Human'

_usage = {
    'AAA': 3.35,
    'AAC': 1.76,
    'AAG': 3.74,
    'AAT': 2.05,
    'ACA': 1.64,
    'ACC': 1.53,
    'ACG': 0.47,
    'ACT': 1.42,
    'AGA': 1.46,
    'AGC': 1.39,
    'AGG': 1.04,
    'AGT': 1.18,
    'ATA': 0.92,
    'ATC': 2.06,
    'ATG': 2.22,
    'ATT': 2.16,
    'CAA': 1.30,
    'CAC': 1.15,
    'CAG': 3.03,
    'CAT': 1.17,
    'CCA': 1.63,
    'CCC': 1.35,
    'CCG': 0.47,
    'CCT': 1.66,
    'CGA': 0.75,
    'CGC': 0.79,
    'CGG': 0.99,
    'CGT': 0.53,
    'CTA': 0.83,
    'CTC': 1.55,
    'CTG': 3.32,
    'CTT': 1.53,
    'GAA': 3.77,
    'GAC': 2.35,
    'GAG': 3.51,
    'GAT': 2.93,
    'GCA': 1.89,
    'GCC': 2.35,
    'GCG': 0.63,
    'GCT': 2.21,
    'GGA': 2.05,
    'GGC': 1.95,
    'GGG': 1.28,
    'GGT': 1.32,
    'GTA': 0.99,
    'GTC': 1.33,
    'GTG': 2.66,
    'GTT': 1.49,
    'TAA': 0.00,
    'TAC': 1.48,
    'TAG': 0.00,
    'TAT': 1.61,
    'TCA': 1.12,
    'TCC': 1.30,
    'TCG': 0.33,
    'TCT': 1.43,
    'TGA': 0.00,
    'TGC': 0.92,
    'TGG': 1.24,
    'TGT': 1.03,
    'TTA': 0.98,
    'TTC': 1.85,
    'TTG': 1.50,
    'TTT': 2.16,
}


class Codon(object):
    """
    A codon consists of three DNA bases, and may be associated with an amino
    # acid. (All codons consisting exclusively of bases ACGT map to an amino or
    stop.)
    """

    def __init__(self, bases):
        """
        :param bases: Any length-three string of IUPAC degenerate bases
            (case-insensitive)
        """
        assert len(bases) == 3, 'codons must have length 3!'
        bases = bases.upper()
        for base in bases:
            assert base in base_utils.ALL_BASES, \
                'unrecognized base "{}"'.format(base)
        self.bases = bases

    def __eq__(self, other):
        """Returns true if two codons have identical bases in the same order."""
        return self.bases == other.bases

    def get_amino(self):
        """
        Returns the string representation of an amino. Typically a
        single-character, except for stop codons. Can return a placeholder for
        an unknown amino in ambiguous cases. (See _UNKNOWN_AMINO)

        :return: string representation of the encoded amino acid.
        """
        options = [base_utils.get_primitives(b) for b in self.bases]
        amino = None
        for base_list in itertools.product(*options):
            base_str = ''.join(base_list)
            if base_str not in _encodings.keys():
                return _UNKNOWN_AMINO
            elif amino and amino != _encodings[base_str]:
                # we've found a conflict
                return _UNKNOWN_AMINO
            else:
                amino = _encodings[base_str]
        return amino

    def encodes_same_amino(self, other):
        """
        Whether two codons encode for the same amino acid.

        Note: '_' is *not* a wildcard. A codon containing this placeholder is
        only considered to amino-match identical codons.
        (e.g. __C encodes the 'same' amino as __C, but not __T, even though this
        substitution would yield the same amino for any specific prefix.)

        :param other: another codon
        :return: boolean
        """
        this_amino = self.get_amino()
        if this_amino == _UNKNOWN_AMINO:
            # unmapped codons must match exactly
            return self.bases == other.bases

        return self.get_amino() == other.get_amino()

    def get_usage(self):
        """
        :return: the usage rate of this codon (genome-dependent)
        """
        return _usage[self.bases] if self.bases in _usage.keys() else 0
