from restriction_enzymes import RestrictionEnzyme
from sequence import Sequence
import restriction_enzymes
import unittest


class TestRestrictionEnzyme(unittest.TestCase):
    def test_init(self):
        enzyme = RestrictionEnzyme('aName', 'CATG')
        self.assertEqual(enzyme.name, 'aName')
        self.assertEqual(enzyme.sequence, Sequence('CATG'))

    def test_len(self):
        self.assertEqual(len(RestrictionEnzyme('x', 'CATG')), 4)
        self.assertEqual(len(RestrictionEnzyme('x', 'GCATGC')), 6)

    def test_is_symmetric(self):
        self.assertTrue(RestrictionEnzyme('x', 'CATG').is_symmetric())
        self.assertFalse(RestrictionEnzyme('x', 'CTTG').is_symmetric())


class TestAllEnzymes(unittest.TestCase):
    def test_returns_a_copy(self):
        first = restriction_enzymes.get_all_enzymes()
        second = restriction_enzymes.get_all_enzymes()
        self.assertEqual(first, second)
        self.assertTrue(first is not second)


if __name__ == '__main__':
    unittest.main()
