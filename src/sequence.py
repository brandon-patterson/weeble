from codons import Codon


class Sequence(object):
    """A sequence of nucleic acid bases. May be unaligned."""

    def __init__(self, base_sequence_str):
        """
        :param base_sequence_str: any case-insensitive string of ACGT_
        """

        # Convert to upper
        base_sequence_str = base_sequence_str.upper()

        # Check that bases are valid
        for base in base_sequence_str:
            assert base in 'ACGT_', 'unrecognized base "{}"'.format(base)

        self.bases = base_sequence_str

    def __eq__(self, other):
        """Whether two Sequences are identical"""
        return self.bases == other.bases

    def __len__(self):
        """The number of bases in the Sequence"""
        return len(self.bases)

    def align(self):
        """
        :return: AlignedSequence obtained by right-padding the Sequence with _'s
        """
        bases_copy = self.bases
        while len(bases_copy) % 3 != 0:
            bases_copy += '_'
        return AlignedSequence(bases_copy)

    def reverse_complement(self):
        """
        Inverts a DNA sequence by flipping each base and reversing the sequence.
        Example: CAT -> ATG

        :returns: the reverse complement Sequence
        """
        return Sequence(self.bases
                        .replace('A', 't')
                        .replace('T', 'a')
                        .replace('C', 'g')
                        .replace('G', 'c')
                        .upper()
                        [::-1])


class AlignedSequence(Sequence):
    """
    A Sequence composed of complete, aligned codons.
    Can be translated into a chain of amino acids.
    """

    def __init__(self, base_sequence):
        """
        :param base_sequence: case-insensitive string of ACGT_, with length 3n
        """
        super().__init__(base_sequence)

        # assert that aligned sequence has proper length
        assert len(base_sequence) % 3 == 0, \
            'AlignedSequence with len {}'.format(len(base_sequence))

        self.codons = []

        for i in range(0, len(base_sequence), 3):
            codon = Codon(base_sequence[i:i + 3])
            self.codons.append(codon)

    def get_amino_string(self):
        """
        :return: string representation of the AlignedSequence's amino chain
        """
        return ''.join(list(map(lambda c: c.get_amino(), self.codons)))

    def encodes_same_aminos(self, other):
        """
        :param other: AlignedSequence
        :return: True iff two sequences represent an identical amino acid chain
        """
        if len(self.codons) != len(other.codons):
            return False
        for left, right in zip(self.codons, other.codons):
            if not left.encodes_same_amino(right):
                return False
        return True
