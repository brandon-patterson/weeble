from sequence import AlignedSequence


class SequenceReplacementEdit(object):
    def __init__(self, aligned_sequence, override_sequence, offset):
        assert offset + len(override_sequence) <= len(aligned_sequence), \
                'override_sequence doesn\'t fit in aligned_sequence with offset'
        self._original_sequence = aligned_sequence
        edited_bases = aligned_sequence.bases[:offset] \
                + override_sequence.bases \
                + aligned_sequence.bases[offset+len(override_sequence):]
        self._new_sequence = AlignedSequence(edited_bases)
        self._edit_begin = offset
        # align to codon begin
        while self._edit_begin % 3 != 0:
            self._edit_begin -= 1
        self._edit_end = offset + len(override_sequence)
        # align to codon end
        while self._edit_end % 3 != 0:
            self._edit_end += 1

    def __eq__(self, other):
        return self._original_sequence == other._original_sequence \
                and self._new_sequence == other._new_sequence

    def __str__(self):
        old = self._original_sequence.bases[self._edit_begin:self._edit_end]
        new = self._new_sequence.bases[self._edit_begin:self._edit_end]
        edit_str = old + ' -> '
        for i in range(len(old)):
            edit_str += new[i] if old[i] == new[i] else new[i].lower()

        return 'edit index: {}\t{}'.format(self._edit_begin+1, edit_str)

    def number_of_bases_modified(self):
        old = self._original_sequence.bases[self._edit_begin:self._edit_end]
        new = self._new_sequence.bases[self._edit_begin:self._edit_end]
        return sum(old[i] != new[i] for i in range(len(old)))

    def number_of_aminos_modified(self):
        begin = self._edit_begin // 3
        end = self._edit_end // 3
        old = self._original_sequence.codons[begin:end]
        new = self._new_sequence.codons[begin:end]
        return sum(not old[i].encodes_same_amino(new[i])
                   for i in range(len(old)))
