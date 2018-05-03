import unittest
import subprocess
import goenrich
from goenrich.ontol import GOntology
import networkx

class TestGoenrich(unittest.TestCase):

    def test_analysis_and_export(self):
        O = GOntology.from_obo('db/go-basic.obo')
        gene2go = goenrich.read.gene2go('db/gene2go.gz')
        values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
        background_attribute = 'gene2go'
        O.propagate(values, background_attribute)
        query = gene2go['GeneID'].unique()[:20]
        try:
            import pygraphviz
            df = O.get_enrichment(query, background_attribute, gvfile='test.dot')
            subprocess.check_call(['dot', '-Tpng', 'test.dot', '-o', 'test.png'])
            subprocess.check_call(['rm', 'test.dot', 'test.png'])
        except ImportError:
            df = O.get_enrichment(query, background_attribute)
            print('pygraphviz could not be imported')
        self.assertEqual(len(df.query('q<0.05')), 8)


    def test_pval_correctness(self):
        G = networkx.DiGraph()
        G.add_node(0, name='r', namespace='a')
        G.add_node(1, name='1', namespace='a')
        G.add_node(2, name='2', namespace='a')
        G.add_edge(0, 1)
        G.add_edge(0, 2)
        O = GOntology.from_graph(G)
        values = {1: set(range(10)), 2: set(range(20))}
        background_attribute = 'bg_attr'
        O.propagate(values, background_attribute)
        query = [1, 2, 3, 4, 5]
        df = O.get_enrichment(query, background_attribute)
        best_pval = float(df.dropna().sort_values('p').head(1).p)
        self.assertAlmostEqual(best_pval, 0.016253869969040255)

if __name__ == '__main__':
    unittest.main()
