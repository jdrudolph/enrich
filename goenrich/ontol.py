from collections import defaultdict
import networkx as nx
from goenrich.obo import ontology as read_obo
import os

class GOTermNotFoundError(Exception):
    def __init__(self, termid, *args, **kwargs):
        self.termid = termid
        super().__init__(*args, **kwargs)

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
    Can be instatiated using the ontology in obo file format.

    >>> O = GOntology('/path/to/obo/file')

    Ontology methods accept GOterm objects as well as ter IDs of the
    form GO:XXXXXXX.  When terms are returned only GOterm objects are
    returned.

    """
    
    def __init__(self, obofile):
        self.graph = read_obo(obofile)
        self.obofile = os.path.abspath(obofile)

    def __repr__(self):
        return '<GO Ontology f:{:s}>'.format(self.obofile)

    def get_term(self, termid):
        """
        Get GOterm object by its id of the form GO:XXXXXXX.
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
        more specific one.  t1 and t2 don't have to be immediate child
        parent pair but can be separated with other term in the GO
        graph.
        """
        self._check_terms(t1, t2)
        t1 = getattr(t1, 'id', t1)
        t2 = getattr(t2, 'id', t2)
        ancestors = nx.descendants(self.graph, t1)
        if (t2 in ancestors):
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
        Search the graph for terms containing a patter in their name
        and return a list of those terms.
        """
        
        ids_found = filter(lambda tid: pattern in self.graph.node[tid]['name'],
                           self.graph.nodes())
        terms_found = [self.get_term(i) for i in sorted(ids_found)]
        return terms_found
