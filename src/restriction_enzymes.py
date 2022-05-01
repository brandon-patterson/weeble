from sequence import Sequence
import base_utils


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


_enzymes = [
    RestrictionEnzyme('AciI', 'CCGC'),
    RestrictionEnzyme('AluI', 'AGCT'),
    RestrictionEnzyme('BamHI', 'GGATCC'),
    RestrictionEnzyme('BsiEI', 'CGRYCG'),
    RestrictionEnzyme('EcoRI', 'GAATTC'),
    RestrictionEnzyme('EcoRV', 'GATATC'),
    RestrictionEnzyme('Fnu4HI', 'GCNGC'),
    RestrictionEnzyme('HaeIII', 'GGCC'),
    RestrictionEnzyme('HgaI', 'GACGC'),
    RestrictionEnzyme('HindIII', 'AAGCTT'),
    RestrictionEnzyme('Hpy99I', 'CGWCG'),
    RestrictionEnzyme('KpnI', 'GGTACC'),
    RestrictionEnzyme('MspA1I', 'CMGCKG'),
    RestrictionEnzyme('NotI', 'GCGGCCGC'),
    RestrictionEnzyme('PstI', 'CTGCAG'),
    RestrictionEnzyme('PvuII', 'CAGCTG'),
    RestrictionEnzyme('SacI', 'GAGCTC'),
    RestrictionEnzyme('SalI', 'GTCGAC'),
    RestrictionEnzyme('Sau3AI', 'GATC'),
    RestrictionEnzyme('ScaI', 'AGTACT'),
    RestrictionEnzyme('SgrAI', 'CRCCGGYG'),
    RestrictionEnzyme('SmaI', 'CCCGGG'),
    RestrictionEnzyme('SpeI', 'ACTAGT'),
    RestrictionEnzyme('SphI', 'GCATGC'),
    RestrictionEnzyme('StuI', 'AGGCCT'),
    RestrictionEnzyme('TaqI', 'TCGA'),
    RestrictionEnzyme('XbaI', 'TCTAGA'),
]


def get_all_enzymes():
    """
    :return: a list of 'default' configured restriction enzymes
    """
    return _enzymes.copy()
