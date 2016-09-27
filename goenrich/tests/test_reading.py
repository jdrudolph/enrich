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
        background = goenrich.read.goa('db/gene_association.goa_ref_human.gz')

    def test_goa_from_file_obj(self):
        with gzip.GzipFile('db/gene_association.goa_ref_human.gz') as f:
            background = goenrich.read.goa(f)
            self.assertFalse(f.closed)

    def test_gene2go(self):
        background = goenrich.read.gene2go('db/gene2go.gz')

    def test_gene2go_from_file_obj(self):
        with gzip.GzipFile('db/gene2go.gz') as f:
            background = goenrich.read.gene2go(f)
            self.assertFalse(f.closed)

if __name__ == '__main__':
    unittest.main()
