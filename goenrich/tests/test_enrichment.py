import unittest
import subprocess
import goenrich
import networkx as nx
nx_version = list(map(int, nx.__version__.split('.')))
if nx_version[0] >= 2 and nx_version[1] >= 4:
    nx.Graph.node = nx.Graph.nodes
networkx = nx

class TestGoenrich(unittest.TestCase):

    def test_analysis_and_export(self):
        O = goenrich.obo.ontology('db/go-basic.obo')
        gene2go = goenrich.read.gene2go('db/gene2go.gz')
        values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
        background_attribute = 'gene2go'
        goenrich.enrich.propagate(O, values, background_attribute)
        query = gene2go['GeneID'].unique()[:20]
        try:
            import pygraphviz
            df = goenrich.enrich.analyze(O, query, background_attribute, gvfile='test.dot')
            subprocess.check_call(['dot', '-Tpng', 'test.dot', '-o', 'test.png'])
            subprocess.check_call(['rm', 'test.dot', 'test.png'])
        except ImportError:
            df = goenrich.enrich.analyze(O, query, background_attribute)
            print('pygraphviz could not be imported')
        self.assertEqual(len(df.query('q<0.05')), 8)

    def test_analysis_with_permutation_fdr(self):
        O = goenrich.obo.ontology('db/go-basic.obo')
        gene2go = goenrich.read.gene2go('db/gene2go.gz')
        values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
        background_attribute = 'gene2go'
        goenrich.enrich.propagate(O, values, background_attribute)
        query = gene2go['GeneID'].unique()[:20]
        df = goenrich.enrich.analyze(O, query, background_attribute, method='permutation',
                permutations = 10)
        # TODO test permuation based FDR calculation
        self.skipTest('NaN handling not finished')
        self.assertEqual(len(df.query('q<0.05')), 8)


    def test_pval_correctness_fdr(self):
        O = networkx.DiGraph()
        O.add_node(0, name='r', namespace='a')
        O.add_node(1, name='1', namespace='a')
        O.add_node(2, name='2', namespace='a')
        O.add_edge(0, 1)
        O.add_edge(0, 2)
        values = {1: set(range(10)), 2: set(range(20))}
        background_attribute = 'bg_attr'
        goenrich.enrich.propagate(O, values, background_attribute)
        query = [1, 2, 3, 4, 5]
        df = goenrich.enrich.analyze(O, query, background_attribute)
        best_pval = float(df.dropna().sort_values('p').head(1).p)
        best_qval = float(df.dropna().sort_values('q').head(1).q)
        self.assertAlmostEqual(best_pval, 0.016253869969040255)
        self.assertAlmostEqual(best_qval, 2 * best_pval)

    def test_pval_correctness_bonferroni(self):
        O = networkx.DiGraph()
        O.add_node(0, name='r', namespace='a')
        O.add_node(1, name='1', namespace='a')
        O.add_node(2, name='2', namespace='a')
        O.add_edge(0, 1)
        O.add_edge(0, 2)
        values = {1: set(range(10)), 2: set(range(20))}
        background_attribute = 'bg_attr'
        goenrich.enrich.propagate(O, values, background_attribute)
        query = [1, 2, 3, 4, 5]
        df = goenrich.enrich.analyze(O, query, background_attribute, method='bonferroni')
        best_pval = float(df.dropna().sort_values('p').head(1).p)
        best_qval = float(df.dropna().sort_values('q').head(1).q)
        self.assertAlmostEqual(best_pval, 0.016253869969040255)
        self.assertAlmostEqual(best_qval, 0.0325077399)

if __name__ == '__main__':
    unittest.main()
