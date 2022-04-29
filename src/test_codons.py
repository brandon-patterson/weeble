from codons import Codon
import codons
import unittest


class TestCodon(unittest.TestCase):
    def test_init_enforces_length(self):
        self.assertRaises(AssertionError, lambda: Codon(''))
        self.assertRaises(AssertionError, lambda: Codon('A'))
        self.assertRaises(AssertionError, lambda: Codon('AA'))
        self.assertRaises(AssertionError, lambda: Codon('AAAA'))
        self.assertRaises(AssertionError, lambda: Codon('AAAAA'))

    def test_init_is_case_insensitive(self):
        self.assertEqual(Codon('AGT'), Codon('agt'))

    def test_init_validates_bases(self):
        # These are all valid
        Codon('AAA')
        Codon('CCC')
        Codon('GGG')
        Codon('TTT')
        Codon('___')  # unknown/placeholder

        # These are all invalid
        self.assertRaises(AssertionError, lambda: Codon('UUU'))
        self.assertRaises(AssertionError, lambda: Codon('BBB'))
        self.assertRaises(AssertionError, lambda: Codon('XXX'))
        self.assertRaises(AssertionError, lambda: Codon('***'))
        self.assertRaises(AssertionError, lambda: Codon('5\'AAA'))
        self.assertRaises(AssertionError, lambda: Codon('3\'AAA'))
        self.assertRaises(AssertionError, lambda: Codon('UUU'))

    def test_eq(self):
        self.assertEqual(Codon('AAA'), Codon('AAA'))
        self.assertNotEqual(Codon('AAA'), Codon('AAG'))
        self.assertEqual(Codon('__T'), Codon('__T'))
        self.assertNotEqual(Codon('AAT'), Codon('__C'))
        self.assertNotEqual(Codon('__T'), Codon('__C'))

    def test_get_amino(self):
        self.assertEqual(Codon('AAA').get_amino(), 'K')
        self.assertEqual(Codon('CAT').get_amino(), 'H')
        self.assertEqual(Codon('TAA').get_amino(), 'Ochre')
        self.assertEqual(Codon('___').get_amino(), codons._UNKNOWN_AMINO)
        self.assertEqual(Codon('_AA').get_amino(), codons._UNKNOWN_AMINO)
        self.assertEqual(Codon('AC_').get_amino(), codons._UNKNOWN_AMINO)

    def test_encodes_same_amino(self):
        self.assertTrue(Codon('AAA').encodes_same_amino(Codon('AAG')))
        self.assertTrue(Codon('CTA').encodes_same_amino(Codon('TTG')))

        # We consider stop codons distinct (optional)
        self.assertFalse(Codon('TAA').encodes_same_amino(Codon('TAG')))

        # '_' is considered a 'missing' base, and maps to an 'UNKNOWN' codon
        self.assertTrue(Codon('TC_').encodes_same_amino(Codon('TC_')))
        self.assertTrue(Codon('__A').encodes_same_amino(Codon('__A')))
        # Note that '_' is *not* a wildcard. (TCN all map to the same amino)
        self.assertFalse(Codon('TC_').encodes_same_amino(Codon('TCA')))


if __name__ == '__main__':
    unittest.main()
