from collections import defaultdict
import networkx as nx
from goenrich.obo import ontology as read_obo
from goenrich.enrich import analyze as _analyze
import os

class GOTermNotFoundError(Exception):
    def __init__(self, termid, *args, **kwargs):
        self.termid = termid
        super().__init__(*args, **kwargs)


def _propagate(O, values, attribute):
    """ Propagate values trough the hierarchy"""
    for n in nx.topological_sort(O):
        current = O.node[n].setdefault(attribute, set())
        current.update(values.get(n, set()))
        for p in O[n]:
            O.node[p].setdefault(attribute, set()).update(current)

class GOterm(object):
    """
    Class implementing GO term.
    """

    def __init__(self, id, name='', definition='', n=None):
        self.id = id
        self.name = name
        self.definition = definition
        self.node = n

    def __repr__(self):
        return '<'+self.id + ' ' + self.name+'>'

    def __str__(self):
        return '<GO Term ' + self.id + '>'

    def __eq__(self, other):
        if ((self.id == other.id) and \
            (self.name == other.name) and \
            (self.definition == other.definition)):
            return True
        else:
            return False

class GOntology(object):
    """
    Class representing the Ontology.
    Can be instatiated from an obo file:

    >>> O = GOntology.from_obo('/path/to/obo/file')

    Or directly from networkx graph:

    >>> O = GOntology.from_graph(G)

    Ontology object methods accepting term as an argument accept
    :class:`GOterm` objects as well as term IDs of the form
    GO\:XXXXXXX. Ontology object methods returning terms return
    objects of :class:`GOterm` class.

    """

    def __init__(self, G, obofile=None):
        self.graph = G
        if obofile:
            self.obofile = os.path.abspath(obofile)
        else:
            self.obofile = None

    @classmethod
    def from_obo(cls, obofile):
        """
        Reads GOntology from obo file privided.
        """
        return cls(read_obo(obofile), obofile)

    @classmethod
    def from_graph(cls, G):
        """
        Creates GOntology instance from a networkx DiGraph.
        """
        return cls(G)

    def __repr__(self):
        if self.obofile:
            return '<GO Ontology f:"{:s}" at {:s}>'.format(self.obofile, hex(id(self)))
        else:
            return '<GO Ontology at {:s}>'.format(hex(id(self)))

    def get_term(self, termid):
        """
        Get :class:`GOterm` object by its id of the form GO\:XXXXXXX.
        """
        if (not isinstance(termid, str)):
            raise TypeError('termid parameter should be a string of the form GO:XXXXXXX')
        if not self.graph.has_node(termid):
            raise GOTermNotFoundError(termid, 'Term with id {:s} not found in ontology {:s}'\
                                      .format(termid, repr(self)))
        n = self.graph.node[termid]
        term = GOterm(id=termid, name=n['name'], n=n)
        return term

    def _check_terms(self, *terms):
        for term in terms:
            if (type(term) == str):
                pass
            elif (isinstance(term, GOterm)):
                term = term.id
            if not self.graph.has_node(term):
                raise GOTermNotFoundError(term, 'Term with id {:s} not found in ontology {:s}'\
                                      .format(term, repr(self)))
        return

    def is_child_parent(self, t1, t2):
        """
        Return True if terms t1 and t2 are related by child parent GO
        relashionship. Return False otherwise.

        Parent term is the less specific one and child term is the
        more specific one.  t1 and t2 don't have to be immediate
        child-parent pair but can be separated by a term chain in the
        GO graph.
        """
        self._check_terms(t1, t2)
        t1 = getattr(t1, 'id', t1)
        t2 = getattr(t2, 'id', t2)
        descendants = nx.descendants(self.graph, t1)
        if (t2 in descendants):
            return True
        else:
            return False

    def get_child_terms(self, term):
        """
        Get immediate child terms of a term
        """

        self._check_terms(term)
        term = getattr(term, 'id', term)
        childids = sorted(self.graph.predecessors(term))
        return [self.get_term(childid) for childid in childids]

    def get_parent_terms(self, term):
        """
        Get immediate parent terms of a term
        """

        self._check_terms(term)
        term = getattr(term, 'id', term)
        parentids = sorted(self.graph.successors(term))
        return [self.get_term(parentid) for parentid in parentids]

    def search_terms_by_name(self, pattern):
        """
        Search the ontology for terms by pattern

        :param pattern: a string to scan for in GO terms' names
        """

        ids_found = filter(lambda tid: pattern in self.graph.node[tid]['name'],
                           self.graph.nodes())
        terms_found = [self.get_term(i) for i in sorted(ids_found)]
        return terms_found

    def propagate(self, values, attribute):
        """ Propagate values trough the hierarchy

        >>> from goenrich.ontol import GOntology
        >>> O = GOntology.from_obo('db/go-basic.obo')
        >>> gene2go = goenrich.read.gene2go('db/gene2go.gz')
        >>> values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
        >>> O.propagate(values, 'gene2go')

        Uses topological sorting of the vertices. Since degrees are
        usually low performance is almost linear time.

        :param values: mapping of nodes to set of ids
        :param attribute: name of the attribute
        """
        _propagate(self.graph, values, attribute)

    def get_enrichment(self, query, background_attr, **kwargs):
        """ run enrichment analysis for query

        >>> from goenrich.ontol import GOntology
        >>> O = GOntology.from_obo('db/go-basic.obo')
        >>> gene2go = goenrich.read.gene2go('db/gene2go.gz')
        >>> values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
        >>> O.propagate(values, 'gene2go')
        >>> df = O.get_enrichment(query, 'gene2go')

        :param query: array like of ids
        :param background_attr: string idicating background attribute to use
        :returns: pandas.DataFrame with results
        """

        return _analyze(self.graph, query, background_attr, **kwargs)
