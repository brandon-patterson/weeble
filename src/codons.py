import base_utils
import itertools

_UNKNOWN_AMINO = '?'

_encodings = None
_usage_source = None
_usage_table = None


def load_encodings(encoding_file='encodings'):
    global _encodings
    if _encodings:
        return
    _encodings = {}
    path = '../configs/{}.txt'.format(encoding_file)
    with open(path, 'r') as file:
        for line in file:
            k, v = line.split(': ')
            _encodings[k] = v.strip()


def load_usage(usage_file_name='human'):
    global _usage_source
    global _usage_table
    if _usage_table:
        # Don't allow the table to be altered post-load
        raise RuntimeError('Already loaded usage table {}'
                           .format(_usage_source))
    _usage_source = usage_file_name
    _usage_table = {}
    path = '../configs/usage_tables/{}.txt'.format(usage_file_name)
    with open(path, 'r') as file:
        for line in file:
            k, v = line.split(': ')
            _usage_table[k] = float(v)


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
        if not _encodings:
            load_encodings()
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
        if not _usage_table:
            load_usage()
        if self.bases in _usage_table.keys():
            return _usage_table[self.bases]
        else:
            return 0
