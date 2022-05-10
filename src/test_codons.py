from codons import Codon
import codons
import itertools
import os
import unittest


def get_usage_options():
    options = set()
    for item in os.listdir('../configs/usage_tables'):
        if item.endswith('.txt'):
            options.add(item.split('.')[0])
    return options


def get_expected_keys():
    b = ['A', 'C', 'G', 'T']
    return set([''.join(p) for p in itertools.product(b, b, b)])


class TestConfigs(unittest.TestCase):

    def test_encoding_keys(self):
        self.assertEqual(set(codons._encodings.keys()), get_expected_keys())

    def test_usage(self):
        max_diff = .005 * 64
        for usage_option in get_usage_options():
            codons._usage_table = None  # clear to avoid crashes
            codons.load_usage(usage_option)

            self.assertEqual(codons._usage_source, usage_option)
            self.assertEqual(set(codons._usage_table.keys()),
                             get_expected_keys())
            usage_sum = sum(codons._usage_table.values())
            self.assertLess(usage_sum, 100 + max_diff)
            self.assertGreater(usage_sum, 100 - max_diff)

        # restore defaults
        codons._usage_table = None  # clear to avoid crashes
        codons.load_usage()


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
        Codon('ACT')
        Codon('GU_')
        Codon('RYK')
        Codon('MBV')
        Codon('DHG')
        Codon('SWN')

        # These are all invalid
        self.assertRaises(AssertionError, lambda: Codon('EEE'))
        self.assertRaises(AssertionError, lambda: Codon('XXX'))
        self.assertRaises(AssertionError, lambda: Codon('***'))
        self.assertRaises(AssertionError, lambda: Codon('.- '))
        self.assertRaises(AssertionError, lambda: Codon('5\'AAA'))
        self.assertRaises(AssertionError, lambda: Codon('3\'AAA'))

    def test_eq(self):
        self.assertEqual(Codon('AAA'), Codon('AAA'))
        self.assertNotEqual(Codon('AAA'), Codon('AAG'))
        self.assertEqual(Codon('__T'), Codon('__T'))
        self.assertNotEqual(Codon('AAT'), Codon('__C'))
        self.assertNotEqual(Codon('__T'), Codon('__C'))

    def test_get_amino(self):
        # standard codons
        self.assertEqual(Codon('AAA').get_amino(), 'K')
        self.assertEqual(Codon('CAT').get_amino(), 'H')
        self.assertEqual(Codon('TAA').get_amino(), 'Ochre')

        # some degenerate codons encode actual aminos
        self.assertEqual(Codon('AAY').get_amino(), 'N')
        self.assertEqual(Codon('ACN').get_amino(), 'T')
        self.assertEqual(Codon('AGR').get_amino(), 'R')
        self.assertEqual(Codon('MGR').get_amino(), 'R')

        # some codons cannot be mapped with certainty
        self.assertEqual(Codon('___').get_amino(), codons._UNKNOWN_AMINO)
        self.assertEqual(Codon('_AA').get_amino(), codons._UNKNOWN_AMINO)
        self.assertEqual(Codon('AC_').get_amino(), codons._UNKNOWN_AMINO)
        self.assertEqual(Codon('NNN').get_amino(), codons._UNKNOWN_AMINO)
        self.assertEqual(Codon('YTS').get_amino(), codons._UNKNOWN_AMINO)

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

    def test_get_usage(self):
        self.assertEqual(Codon('ACT').get_usage(), 1.42)
        self.assertEqual(Codon('NBH').get_usage(), 0)


if __name__ == '__main__':
    unittest.main()
