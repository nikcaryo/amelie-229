
from collections import defaultdict

import pkg_resources


def parent_package_name():
    return '.'.join(__name__.split('.')[:-2])


def parse_identifiers(resource):
    rv = defaultdict(lambda: set())
    for line in resource.splitlines():
        line = line.rstrip('\n').split('\t')
        assert len(line) == 2, line
        rv[line[0]].add(line[1])
    return dict(rv)


class HPO:

    def __init__(self, hpo_owl=None, xref_synonyms=False):
        if hpo_owl is None:
            hpo_file = pkg_resources.resource_string(parent_package_name(),
                                                     'pheno_data/hpo_ontology.obo').decode('utf8')
        else:
            hpo_file = pkg_resources.resource_string(parent_package_name(),
                                                     hpo_owl).decode('utf8')
        terms = parse_hpo_text(hpo_file)
        children_mapping = defaultdict(set)
        parents_mapping = defaultdict(set)

        for term in terms.values():
            if 'is_a' in term:
                for parent in term['is_a']:
                    parent_code = parent.split(' ')[0]
                    parents_mapping[term['id'][0]].add(parent_code)
                    children_mapping[parent_code].add(term['id'][0])

        self.all_descendants_mapping = {}
        self.all_ancestors_mapping = {}

        def get_all_children(term_id):
            if term_id not in self.all_descendants_mapping:
                all_descendants = set()
                for child in children_mapping[term_id]:
                    all_descendants.update(get_all_children(child))
                    all_descendants.add(term_id)

                self.all_descendants_mapping[term_id] = all_descendants

            return self.all_descendants_mapping[term_id]
            
        get_all_children('HP:0000118')  # DO NOT TAKE HP:0000001

        self.terms = {}
        self.children_mapping = {}

        def get_all_parents(term_id):
            if term_id not in self.all_ancestors_mapping:
                if term_id not in self.terms:
                    return set()
                all_ancestors = set()
                for parent in parents_mapping[term_id]:
                    if parent not in self.terms:
                        continue
                    all_ancestors.update(get_all_parents(parent))
                    all_ancestors.add(term_id)

                self.all_ancestors_mapping[term_id] = all_ancestors

            return self.all_ancestors_mapping[term_id]

        for term in self.all_descendants_mapping:
            self.terms[term] = terms[term]
            self.children_mapping[term] = children_mapping[term]

        for term in self.terms:
            get_all_parents(term)

        if xref_synonyms:
            ontology_mappings = {'UMLS': parse_identifiers(pkg_resources.resource_string(parent_package_name(),
                                              'pheno_data/hpo_umls_terms.tsv').decode('utf8')),
                                 'SNOMEDCT_US': parse_identifiers(
                                     pkg_resources.resource_string(parent_package_name(),
                                                                   'pheno_data/hpo_snomed_terms.tsv').decode('utf8')),
                                 'MSH': parse_identifiers(
                                     pkg_resources.resource_string(parent_package_name(),
                                                                   'pheno_data/hpo_mesh_terms.tsv').decode('utf8'))}

            for term_id in self.terms:
                if 'xref' in self.terms[term_id]:
                    xrefs = self.terms[term_id]['xref']
                    for xref in xrefs:
                        xref = xref.split(':')
                        if len(xref) < 2:
                            continue
                        ontology, identifier = xref[:2]
                        if ontology in ontology_mappings:
                            if 'xref_synonyms' not in self.terms[term_id]:
                                self.terms[term_id]['xref_synonyms'] = set()
                            if identifier in ontology_mappings[ontology]:
                                for synonym in ontology_mappings[ontology][identifier]:
                                    # print("Adding xref synonym %s from %s to HPO term %s ")
                                    self.terms[term_id]['xref_synonyms'].add(synonym)

    def __iter__(self):
       yield from self.terms

    def __len__(self):
        return len(self.terms)

    def get_name(self, term_id):
        return self.terms[term_id]['name'][0]

    def get_synonyms(self, term_id):
        result = []

        result.append(self.terms[term_id]['name'][0])

        if 'synonym' in self.terms[term_id]:
            for synonym in self.terms[term_id]['synonym']:
                string_part = synonym.split('"', 3)[1]
                result.append(string_part)

        return result

    def get_xref_synonyms(self, term_id):
        result = []
        if 'xref_synonyms' in self.terms[term_id]:
            for synonym in self.terms[term_id]['xref_synonyms']:
                result.append(synonym)
        return result

    def has_descendant(self, term_id1, term_id2):
        return term_id2 in self.all_descendants_mapping[term_id1]

    def canonicalize(self, hpo_ids):
        rv = set()
        for hpo_id in hpo_ids:
            if hpo_id not in self.terms:
                continue
            rv.update(self.all_ancestors_mapping[hpo_id])
        return rv - {'HP:0000118'}


def parse_hpo_text(text):
    terms = {}

    current_term = None
    for line in text.splitlines():
        if line == '[Term]':
            if current_term is not None:
                terms[current_term['id'][0]] = dict(current_term)

            current_term = defaultdict(list)
        elif current_term is not None:
            if ':' not in line:
                continue

            key, value = line.split(': ', 1)
            current_term[key].append(value)

    if current_term is not None:
        terms[current_term['id'][0]] = dict(current_term)

    return terms