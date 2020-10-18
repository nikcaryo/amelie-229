import pkg_resources
import io
import re
from collections import defaultdict
from .clinvar import parent_package_name

from amelie import tuples

no_alnum = re.compile(r'[^a-zA-Z0-9_]+')


def normalize_gene_name(symbol, mapping_type):
    symbol = [no_alnum.sub('', a) for a in symbol]
    # this costs a TON of time and doesn't actually get us that much closer, I think ...
    # symbol = replace_some_greek_letters_list(symbol)
    if not mapping_type.startswith('PROTEIN'):
        return ' '.join(symbol)
    if len(symbol) == 1 and any(x.isupper() for x in symbol[0][1:]):
        return symbol[0]
    for j in range(len(symbol)):
        symbol[j] = symbol[j].lower()
    symbol = ' '.join(symbol)
    rv = ' '.join(no_alnum.sub(' ', symbol).split())
    return rv


def gene_symbol_to_ensembl_id_map():
    lines = pkg_resources.resource_string(parent_package_name(), 'gene_data/ensembl_genes_with_names.tsv').\
        decode('utf8').strip().split('\n')
    eid_map = defaultdict(set)
    eid_map_permuted = defaultdict(set)
    for line in lines:
        eid, canonical_name, orig_name, mapping_type = line.rstrip('\n').split('\t')
        orig_name = orig_name.strip()
        gene_name = orig_name.split()
        gene_name = normalize_gene_name(gene_name, mapping_type)
        eid_map[orig_name].add(tuples.GeneName(eid, canonical_name, mapping_type))
        eid_map[gene_name].add(tuples.GeneName(eid, canonical_name, mapping_type))
        eid_map[gene_name.replace(' ', '')].add(tuples.GeneName(eid, canonical_name, mapping_type))
        if mapping_type.startswith('PROTEIN'):
            eid_map_permuted[frozenset(gene_name.split())].add(tuples.GeneName(eid, canonical_name, mapping_type))
    return eid_map, eid_map_permuted


def ensembl_id_map_to_gene_symbol():
    lines = pkg_resources.resource_string(parent_package_name(), 'gene_data/ensembl_genes_with_names.tsv').\
        decode('utf8').strip().split('\n')
    eid2symbol = {}
    for line in lines:
        eid, symbol, modified_symbol, col_type = line.strip().split("\t")
        eid2symbol[eid] = symbol
    return eid2symbol


def load_ensembl_to_entrez():
    ensembl_to_entrez = {}
    with pkg_resources.resource_stream(parent_package_name(), 'gene_data/gene.entrez.map') as mapData:
        with io.TextIOWrapper(mapData) as textMapData:
            for line in textMapData:
                line = line.strip().split()
                ensembl_id = line[0]
                entrez_id = line[1]
                ensembl_to_entrez[ensembl_id] = entrez_id

    return ensembl_to_entrez


def load_entrez_to_ensembl():
    entrez_to_ensembl = {}
    with pkg_resources.resource_stream(parent_package_name(), 'gene_data/gene.entrez.map') as mapData:
        with io.TextIOWrapper(mapData) as textMapData:
            for line in textMapData:
                line = line.strip().split()
                ensembl_id = line[0]
                entrez_id = line[1]
                entrez_to_ensembl[entrez_id] = ensembl_id
    return entrez_to_ensembl
