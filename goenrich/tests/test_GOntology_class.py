import unittest
from pkg_resources import resource_filename as pkg_file
from goenrich.ontol import GOntology, GOTermNotFoundError

class TestGOntology(unittest.TestCase):
    O = GOntology(pkg_file(__name__, 'test_ontologies/goslim_generic.obo'))

    def test_creation(self):
        self.assertTrue('test_ontologies/goslim_generic.obo' in self.O.obofile)
        e = self.O.graph.get_edge_data('GO:0000228', 'GO:0005634')
        self.assertEqual(e['rel'], 'part_of')

    def test_get_term(self):
        t = self.O.get_term('GO:0003729')
        self.assertEqual(t.name, 'mRNA binding')
        with self.assertRaises(GOTermNotFoundError):
            self.O.get_term('GO:1234567')

    def test_is_child_parent(self):
        self.assertTrue(self.O.is_child_parent('GO:0065003', 'GO:0022607'))
        self.assertFalse(self.O.is_child_parent('GO:0065003', 'GO:0000902'))

    def test_get_child_terms(self):
        self.assertListEqual(self.O.get_child_terms('GO:0009056'),
                             list(map(lambda t: self.O.get_term(t),
                                      ['GO:0006914', 'GO:0034655'])))
    def test_get_parent_terms(self):
        self.assertListEqual(self.O.get_parent_terms('GO:0034655'),
                             list(map(lambda t: self.O.get_term(t),
                                      ['GO:0009056', 'GO:0034641'])))

    def test_search_term_by_name(self):
        self.assertListEqual(self.O.search_terms_by_name('mitochond'),
                             list(map(lambda t: self.O.get_term(t),
                                      ['GO:0005739', 'GO:0007005'])))

if __name__ == '__main__':
    unittest.main()
