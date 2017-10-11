import pkg_resources
import unittest
import gzip
import goenrich

class TestRead(unittest.TestCase):
    def test_ontology(self):
        G = goenrich.obo.ontology('db/go-basic.obo')

    def test_ontology_from_file_obj(self):
        with open('db/go-basic.obo') as f:
            G = goenrich.obo.ontology(f)
            self.assertFalse(f.closed)

    def test_goa(self):
        background = goenrich.read.goa('db/gene_association.goa_human.gaf.gz')

    def test_goa_from_file_obj(self):
        with gzip.GzipFile('db/gene_association.goa_human.gaf.gz') as f:
            background = goenrich.read.goa(f)
            self.assertFalse(f.closed)

    def test_gene2go(self):
        background = goenrich.read.gene2go('db/gene2go.gz')

    def test_gene2go_from_file_obj(self):
        with gzip.GzipFile('db/gene2go.gz') as f:
            background = goenrich.read.gene2go(f)
            self.assertFalse(f.closed)

    def test_goslim_from_file(self):
        G = goenrich.obo.ontology(pkg_resources.resource_filename(goenrich.__name__, 'tests/test_ontologies/goslim_generic.obo'))
        self.assertEqual(len(G.nodes()), 150)
        self.assertSetEqual(set(G.successors('GO:0009056')), set(['GO:0008150']))
        self.assertSetEqual(set(G.predecessors('GO:0009056')), set(['GO:0034655', 'GO:0006914']))

if __name__ == '__main__':
    unittest.main()
