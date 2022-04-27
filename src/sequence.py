from codons import Codon


class Sequence(object):
    def __init__(self, base_sequence_str):
        # Convert to upper
        base_sequence_str = base_sequence_str.upper()

        # Check that bases are valid
        for base in base_sequence_str:
            assert base in 'ACGT_', 'unrecognized base "{}"'.format(base)

        self.bases = base_sequence_str

    def __eq__(self, other):
        return self.bases == other.bases

    def __len__(self):
        return len(self.bases)

    def align(self):
        bases_copy = self.bases
        while len(bases_copy) % 3 != 0:
            bases_copy += '_'
        return AlignedSequence(bases_copy)

    def invert(self):
        """
        Inverts a DNA sequence by flipping each base and reversing the sequence
        (CAT -> ATG)
        """
        return Sequence(self.bases
                        .replace('A', 't')
                        .replace('T', 'a')
                        .replace('C', 'g')
                        .replace('G', 'c')
                        .upper()
                        [::-1])


class AlignedSequence(Sequence):
    def __init__(self, base_sequence):
        super().__init__(base_sequence)

        # assert that aligned sequence has proper length
        assert len(base_sequence) % 3 == 0, \
            'AlignedSequence with len {}'.format(len(base_sequence))

        self.codons = []

        for i in range(0, len(base_sequence), 3):
            codon = Codon(base_sequence[i:i+3])
            self.codons.append(codon)

    def get_amino_string(self):
        return ''.join(list(map(lambda c: c.get_amino(), self.codons)))

    def encodes_same_aminos(self, other):
        if len(self.codons) != len(other.codons):
            return False
        for left, right in zip(self.codons, other.codons):
            if not left.encodes_same_amino(right):
                return False
        return True
