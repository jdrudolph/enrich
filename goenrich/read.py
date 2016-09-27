"""
parsers for different go-annotation formats
"""
import pandas as pd
GENE_ASSOCIATION_COLUMNS = ('db', 'db_object_id', 'db_object_symbol',
                            'qualifier', 'go_id', 'db_reference',
                            'evidence_code', 'with_from', 'aspect',
                            'db_object_name', 'db_object_synonym',
                            'db_object_type', 'taxon', 'date', 'assigned_by',
                            'annotation_extension', 'gene_product_form_id')
EXPERIMENTAL_EVIDENCE = ('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP')
def goa(filename, experimental=True, **kwds):
    """ read go-annotation file
    
    :param filename: protein or gene identifier column
    :param experimental: use only experimentally validated annotations
    """
    defaults = {'comment' : '!',
            'names': GENE_ASSOCIATION_COLUMNS}

    if experimental and 'usecols' in kwds:
        kwds['usecols'] += ('evidence_code', )

    defaults.update(kwds)
    result = pd.read_table(filename, **defaults)

    if experimental:
        retain_mask = result.evidence_code.isin(EXPERIMENTAL_EVIDENCE)
        result.drop(result.index[~retain_mask], inplace=True)

    return result

def sgd(filename, experimental=False, **kwds):
    """ read yeast genome database go-annotation file

    :param filename: protein or gene identifier column
    :param experimental: use only experimentally validated annotations
    """
    return goa(filename, experimental, **kwds)

GENE2GO_COLUMNS = ('tax_id', 'GeneID', 'GO_ID', 'Evidence', 'Qualifier', 'GO_term', 'PubMed', 'Category')
def gene2go(filename, experimental=False, tax_id=9606, **kwds):
    """ read go-annotation file
        
    :param filename: protein or gene identifier column
    :param experimental: use only experimentally validated annotations
    :param tax_id: filter according to taxon
    """
    defaults = {'comment': '#',
                'names': GENE2GO_COLUMNS}
    defaults.update(kwds)
    result = pd.read_table(filename, **defaults)
    
    retain_mask = result.tax_id == tax_id
    result.drop(result.index[~retain_mask], inplace=True)

    if experimental:
        retain_mask = result.Evidence.isin(EXPERIMENTAL_EVIDENCE)
        result.drop(result.index[~retain_mask], inplace=True)

    return result
