from sequence_edit import SequenceReplacementEdit


class WobbleCutDetector(object):
    """Detects wobble-enabled restriction enzyme cuts"""

    def __init__(self):
        pass

    def detect_cuts(self, aligned_sequence, restriction_enzyme):
        """
        Returns a list of potential edits that enable the aligned sequence to
        be cut by a given restriction enzyme. Includes cuts matching the
        reverse complement and preexisting cut-sites (e.g. zero-base edits).

        :param aligned_sequence: AlignedSequence
        :param restriction_enzyme: RestrictionEnzyme
        :return: List of SequenceReplacementEdit enabling an enzyme cut
        """
        edit_list = self._detect_cuts_one_way(
            aligned_sequence, restriction_enzyme.sequence)
        if not restriction_enzyme.is_symmetric():
            edit_list += self._detect_cuts_one_way(
                aligned_sequence,
                restriction_enzyme.sequence.reverse_complement())
        return edit_list

    @staticmethod
    def _detect_cuts_one_way(aligned_sequence, cut_seq):
        """
        Helper that detects edit-enabled cut-sites for the given inputs (without
        considering the reverse complement).

        :param aligned_sequence: AlignedSequence
        :param cut_seq: Sequence from a RestrictionEnzyme of interest
        :return: List of SequenceReplacementEdits enabling a cut
        """
        edit_list = []
        for offset in range(len(aligned_sequence) - len(cut_seq) + 1):
            edit = SequenceReplacementEdit(aligned_sequence, cut_seq, offset)
            if edit.number_of_aminos_modified() == 0:
                edit_list.append(edit)
        return edit_list
