import base_utils
import unittest


class TestBaseUtils(unittest.TestCase):
    def test_get_primitives(self):
        self.assertEqual(base_utils.get_primitives('A'), ['A'])
        self.assertEqual(base_utils.get_primitives('C'), ['C'])
        self.assertEqual(base_utils.get_primitives('G'), ['G'])
        self.assertEqual(base_utils.get_primitives('T'), ['T'])
        self.assertEqual(base_utils.get_primitives('_'), ['_'])
        self.assertEqual(base_utils.get_primitives('a'), ['A'])
        self.assertEqual(base_utils.get_primitives('U'), ['T'])
        self.assertEqual(base_utils.get_primitives('R'), ['A', 'G'])
        self.assertEqual(base_utils.get_primitives('B'), ['C', 'G', 'T'])
        self.assertEqual(base_utils.get_primitives('N'), ['A', 'C', 'G', 'T'])


if __name__ == '__main__':
    unittest.main()
