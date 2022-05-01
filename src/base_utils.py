PRIMITIVE_BASES = frozenset(['A', 'C', 'G', 'T', 'U', '_'])
DEGENERATE_BASES = frozenset(['B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'Y'])
ALL_BASES = frozenset.union(PRIMITIVE_BASES, DEGENERATE_BASES)


def get_primitives(base):
    """
    :param base: a single-character base representation (case insensitive)
    :return: a list of all primitive bases the input can represent (excluding U)
    """
    base = base.upper()
    if base in PRIMITIVE_BASES:
        return [base] if base != 'U' else ['T']
    if base == 'K':
        return ['G', 'T']
    if base == 'M':
        return ['A', 'C']
    if base == 'R':
        return ['A', 'G']
    if base == 'S':
        return ['C', 'G']
    if base == 'W':
        return ['A', 'T']
    if base == 'Y':
        return ['C', 'T']
    if base == 'B':
        return ['C', 'G', 'T']
    if base == 'D':
        return ['A', 'G', 'T']
    if base == 'H':
        return ['A', 'C', 'T']
    if base == 'V':
        return ['A', 'C', 'G']
    if base == 'N':
        return ['A', 'C', 'G', 'T']
