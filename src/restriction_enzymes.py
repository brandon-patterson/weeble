from sequence import Sequence


class RestrictionEnzyme(object):
    def __init__(self, name, base_sequence):
        self.name = name
        self.sequence = Sequence(base_sequence)

    def __len__(self):
        return len(self.sequence)

    def is_symmetric(self):
        return self.sequence == self.sequence.invert()


_enzymes = [
    RestrictionEnzyme('AciI', 'CCGC'),
    RestrictionEnzyme('AluI', 'AGCT'),
    RestrictionEnzyme('BamHI', 'GGATCC'),
    RestrictionEnzyme('EcoRI', 'GAATTC'),
    RestrictionEnzyme('EcoRV', 'GATATC'),
    RestrictionEnzyme('HaeIII', 'GGCC'),
    RestrictionEnzyme('HgaI', 'GACGC'),
    RestrictionEnzyme('HindIII', 'AAGCTT'),
    RestrictionEnzyme('KpnI', 'GGTACC'),
    RestrictionEnzyme('NotI', 'GCGGCCGC'),
    RestrictionEnzyme('PstI', 'CTGCAG'),
    RestrictionEnzyme('PvuII', 'CAGCTG'),
    RestrictionEnzyme('SacI', 'GAGCTC'),
    RestrictionEnzyme('SalI', 'GTCGAC'),
    RestrictionEnzyme('Sau3AI', 'GATC'),
    RestrictionEnzyme('ScaI', 'AGTACT'),
    RestrictionEnzyme('SmaI', 'CCCGGG'),
    RestrictionEnzyme('SpeI', 'ACTAGT'),
    RestrictionEnzyme('SphI', 'GCATGC'),
    RestrictionEnzyme('StuI', 'AGGCCT'),
    RestrictionEnzyme('TaqI', 'TCGA'),
    RestrictionEnzyme('XbaI', 'TCTAGA'),
]


def get_all_enzymes():
    return _enzymes.copy()
