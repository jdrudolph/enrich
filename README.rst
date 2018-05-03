goenrich
========

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/jdrudolph/goenrich?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

.. image:: https://readthedocs.org/projects/goenrich/badge/?version=latest
    :target: https://goenrich.readthedocs.org/en/latest

.. image:: https://travis-ci.org/jdrudolph/goenrich.svg?branch=master
    :target: https://travis-ci.org/jdrudolph/goenrich

Convenient GO enrichments from python. For use in ``python`` projects.

#. Builds the GO-ontology graph
#. Propagates GO-annotations up the graph
#. Performs enrichment test for all categories
#. Performs multiple testing correction
#. Allows for export to ``pandas`` for processing and ``graphviz`` for
   visualization

Installation
------------

| Install package from pypi and download ontology and needed annotations.

.. code:: shell

    pip install goenrich
    # assuming in goenrich distribution root directory (containing setup.py file)
    mkdir db
    # Ontology
    wget http://purl.obolibrary.org/obo/go/go-basic.obo -O db/go-basic.obo
    # UniprotACC
    wget http://geneontology.org/gene-associations/goa_human.gaf.gz -O db/gene_association.goa_human.gaf.gz
    # Yeast SGD
    wget http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gz -O db/gene_association.sgd.gz
    # Entrez GeneID
    wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz -O db/gene2go.gz

Usage
-----

Create and use the Ontology
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

  import goenrich
  import goenrich.ontol

  # build the ontology
  O = goenrich.ontol.GOntology.from_obo('db/go-basic.obo')

  # getting a term
  t = O.get_term('GO:0003729')
  # <GO:0003729 mRNA binding>

  # getting GO term child terms
  tt = O.get_child_terms(t) #or with Term ID O.get_child_terms('GO:0003729')
  # [<GO:0003730 mRNA 3'-UTR binding>,
  #  <GO:0030350 iron-responsive element binding> ...]

  # search GO term by name
  O.search_terms_by_name('translation factor activity')
  # [<GO:0008135 translation factor activity, RNA binding>,
  #  <GO:0045183 translation factor activity, non-nucleic acid binding>]

See docs of :py:class:`~goenrich.ontol.GOntology` class for more goodies.

Run GO enrichment
^^^^^^^^^^^^^^^^^

.. code:: python

  # use all entrez geneid associations form gene2go as background
  # use annot = goenrich.read.goa('db/gene_association.goa_human.gaf.gz') for uniprot
  # use annot = goenrich.read.sgd('db/gene_association.sgd.gz') for yeast
  gene2go = goenrich.read.gene2go('db/gene2go.gz')
  # use values = {k: set(v) for k,v in annot.groupby('go_id')['db_object_symbol']} for uniprot/yeast
  values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}

  # propagate the background through the ontology
  background_attribute = 'gene2go'
  O.propagate(values, background_attribute)

  # extract some list of entries as example query
  # use query = annot['db_object_symbol'].unique()[:20]
  query = gene2go['GeneID'].unique()[:20]

  # for additional export to graphviz just specify the gvfile argument
  # the show argument keeps the graph reasonably small
  df = O.get_enrichment(query, background_attribute, gvfile='test.dot')

  # generate html
  df.dropna().head().to_html('example.html')

  # call to graphviz
  import subprocess
  subprocess.check_call(['dot', '-Tpng', 'test.dot', '-o', 'test.png'])

Generate ``png`` image using graphviz:

.. code:: shell

    dot -Tpng example.dot > example.png

or directly from python:

.. code:: python

  import subprocess
  subprocess.check_call(['dot', '-Tpng', 'example.dot', '-o', 'example.png'])

.. image:: https://cloud.githubusercontent.com/assets/2606663/8525018/cad3a288-23fe-11e5-813c-bd205a47eed8.png

Check the documentation for all available parameters

Licence & Contributors
======================

This work is licenced under the MIT licence

Contributions are welcome!

Special thanks
--------------

- `@lukauskas <https://github.com/lukauskas/>`_ for implementing i/o support for file-like objects.
- `@zfrenchee <https://github.com/zfrenchee/>`_ for fixing a bug in the calculation of the test statistic.
- `@mikpom <https://github.com/mikpom/>`_ for implementing support for `networkx >= 2.0.0`.

Building the documentation
==========================

.. code:: shell

  sphinx-apidoc -f -o docs goenrich goenrich/tests
  sphinx-build docs docs/_build
