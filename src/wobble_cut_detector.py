from sequence_edit import SequenceReplacementEdit


class WobbleCutDetector(object):
    def __init__(self):
        pass

    def detect_cuts(self, aligned_sequence, restriction_enzyme):
        edit_list = self._detect_cuts_one_way(
            aligned_sequence, restriction_enzyme.sequence)
        if not restriction_enzyme.is_symmetric():
            edit_list += self._detect_cuts_one_way(
                aligned_sequence, restriction_enzyme.sequence.invert())
        return edit_list

    @staticmethod
    def _detect_cuts_one_way(aligned_sequence, cut_seq):
        edit_list = []
        for offset in range(len(aligned_sequence) - len(cut_seq) + 1):
            edit = SequenceReplacementEdit(aligned_sequence, cut_seq, offset)
            if edit.number_of_aminos_modified() == 0:
                edit_list.append(edit)
        return edit_list
