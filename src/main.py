"""
Author: Brandon Patterson

A self-contained toy script for determining possible cuts of a nucleotide
sequence, including potential wobble-enabled cuts.
"""

# Must be properly aligned. (Use leading underscores to align first codon.)
sequence_of_interest = 'AGTGGCTCCTAGAGGCTCGAATCAGCTATC'

# TODO: update enzymes based on ~availability
# Sorting by preference will list common enzymes first in results
restriction_enzymes = {
    'AciI': 'CCGC',
    'AluI': 'AGCT',
    'BamHI': 'GGATCC',
    'EcoRI': 'GAATTC',
    'EcoRV': 'GATATC',
    'HaeIII': 'GGCC',
    'HgaI': 'GACGC',
    'HindIII': 'AAGCTT',
    'KpnI': 'GGTACC',
    'NotI': 'GCGGCCGC',
    'PstI': 'CTGCAG',
    'PvuII': 'CAGCTG',
    'SacI': 'GAGCTC',
    'SalI': 'GTCGAC',
    'Sau3AI': 'GATC',
    'ScaI': 'AGTACT',
    'SmaI': 'CCCGGG',
    'SpeI': 'ACTAGT',
    'SphI': 'GCATGC',
    'StuI': 'AGGCCT',
    'TaqI': 'TCGA',
    'XbaI': 'TCTAGA',
}

# TODO: double check these, they were copied manually...
codons = {
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


def equivalent_codons(codon_one, codon_two):
    """
    Returns True if two codons are guaranteed to code the same amino.
    ('_' represents an unknown amino)
    """
    if codon_one == codon_two:
        return True
    if '_' in codon_one or '_' in codon_two:
        # TODO: consider optimizing (e.g. __T and __C *are* equivalent)
        return False
    return codons[codon_one] == codons[codon_two]


def equivalent_sequences(seq_one, seq_two):
    """
    Returns True if two codon sequences code for the same amino sequence.
    ('_' represents an unknown amino)
    """
    if len(seq_one) != len(seq_two):
        return False

    assert (len(seq_one) % 3 == 0)

    # for codon in seq
    for i in range(0, len(seq_one), 3):
        if not equivalent_codons(seq_one[i:i + 3], seq_two[i:i + 3]):
            return False

    return True


def invert(seq):
    """
    Inverts a DNA sequence by flipping each base and reversing the sequence
    (CAT -> ATG)
    """
    return seq.replace('A', 't').replace('T', 'a').replace('C', 'g') \
               .replace('G', 'c').upper()[::-1]


def get_edit(seq_from, seq_to):
    """
    Returns the 'edit' sequence needed to convert seq_from to seq_to.
    For example, get_edit('ATGGTC', 'ATCGAC') returns 'ATcGaC'
    """
    assert (len(seq_from) == len(seq_to))

    edit = ''
    for i in range(len(seq_from)):
        if seq_from[i] == seq_to[i]:
            edit += seq_to[i]
        else:
            edit += seq_to[i].lower()

    return edit


def detect_wobble_cuts(seq, cut_seq):
    """
    Returns a list of location/edit pairs that will enable a cut_seq match.
    An 'edit' is a copy of the cut sequence with altered aminos in lowercase.
    """
    edit_list = []
    cut_set = set()
    cut_set.add(cut_seq)
    cut_set.add(invert(cut_seq))
    for pattern in cut_set:
        for i in range(len(seq) - len(pattern) + 1):
            seq_edit = seq[:i] + pattern + seq[i + len(pattern):]
            # TODO: compare relevant portions of sequence to improve performance
            if equivalent_sequences(seq, seq_edit):
                edit = get_edit(seq[i:i + len(pattern)], pattern)
                edit_list.append((i, edit))
    return edit_list


if __name__ == '__main__':
    print('Checking sequence for potential cuts...')
    print(sequence_of_interest)
    print()

    # pad sequence to a full frame
    padded_sequence = sequence_of_interest
    while len(padded_sequence) % 3 != 0:
        padded_sequence += '_'

    for name, cut_pattern in restriction_enzymes.items():
        wobble_cuts = detect_wobble_cuts(padded_sequence, cut_pattern)
        if len(wobble_cuts) > 0:
            print('possible cuts for {}:'.format(name))
            for cut in wobble_cuts:
                print('  index: {}  edit: {}'.format(cut[0], cut[1]))

            print()
