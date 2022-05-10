from sequence import Sequence

_enzymes = None


class RestrictionEnzyme(object):
    """
    A restriction enzyme, which will cut any DNA sequence matching a pattern.
    """

    def __init__(self, name, base_sequence):
        """
        :param name: the name of the enzyme
        :param base_sequence: a string of bases (of any length)
        """
        self.name = name
        self.sequence = Sequence(base_sequence)

    def __len__(self):
        """The number of bases in the match pattern."""
        return len(self.sequence)

    def is_symmetric(self):
        """
        Whether the enzyme's match pattern is its own reverse complement.

        Examples:
            AATT -> True
            GGTACC -> True
            GCCG -> False

        :return: boolean
        """
        return self.sequence == self.sequence.reverse_complement()


def load_enzymes(enzyme_file):
    global _enzymes
    if _enzymes:
        # Prevent reloading of enzymes mid-run, which could cause errors.
        raise RuntimeError('enzymes already loaded')
    _enzymes = []
    with open('../configs/{}.txt'.format(enzyme_file), 'r') as infile:
        for line in infile:
            name, pattern = line.split(': ')
            _enzymes.append(RestrictionEnzyme(name, pattern.strip()))


def get_all_enzymes():
    """
    :return: a list of 'default' configured restriction enzymes
    """
    if not _enzymes:
        load_enzymes('restriction_enzymes')
    return _enzymes.copy()
