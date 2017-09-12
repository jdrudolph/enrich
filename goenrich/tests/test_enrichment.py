import unittest
import subprocess
import goenrich
import networkx

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
            goenrich.enrich.analyze(O, query, background_attribute, gvfile='test.dot')
            subprocess.check_call(['dot', '-Tpng', 'test.dot', '-o', 'test.png'])
            subprocess.check_call(['rm', 'test.dot', 'test.png'])
        except ImportError:
            goenrich.enrich.analyze(O, query, background_attribute)
            print('pygraphviz could not be imported')

    def test_pval_correctness(self):
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
        print(best_pval)
        assert best_pval == float(0.016253869969040255)

if __name__ == '__main__':
    unittest.main()
