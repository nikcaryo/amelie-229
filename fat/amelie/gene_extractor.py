import re
from collections import defaultdict
import unicodedata

from . import tuples
from amelie.dataloaders.gene_data import normalize_gene_name, gene_symbol_to_ensembl_id_map

import traceback


ensembl_mapping_types = ['HGNC_SYMBOL',
                         'PROTEIN_GENE2PUBTATOR',
                         'HGNC_SYNONYM',
                         'MANUAL_SYNONYM',
                         'PROTEIN_FULL_NAME',
                         'PROTEIN_SHORT_NAME',
                         'UNIPROT_GENE_NAME', ]

max_len = 8

min_word_len = {
  'HGNC_SYMBOL': 2,
  'HGNC_SYNONYM': 3,
  'MANUAL_SYNONYM': 3,
  'PROTEIN_FULL_NAME': 5,
  'PROTEIN_SHORT_NAME': 4,
  'UNIPROT_GENE_NAME': 3,
  'PROTEIN_GENE2PUBTATOR': 3
}


def rgx_mult_search(phrase, strings, rgxs, orig_strings, orig_rgxs, flags=re.I):
    for i, s in enumerate(strings):
        try:
            regex = re.escape(s)
            if re.search(regex, phrase, flags):
                return orig_strings[i]
        except Exception:
            traceback.print_exc()
            print(regex)
    for i, s in enumerate(rgxs):
        try:
            regex = s
            if re.search(regex, phrase, flags):
                return orig_rgxs[i]
        except Exception:
            traceback.print_exc()
            print(regex)
    return None


def replace_some_greek_letters(s):
    s = s.replace('alpha', 'a')
    s = s.replace('beta', 'b')
    s = s.replace('gamma', 'g')
    s = s.replace('α', 'a')
    s = s.replace('𝛼', 'a')
    s = s.replace('β', 'b')
    s = s.replace('γ', 'g')
    return ''.join(c for c in unicodedata.normalize('NFD', s)
                   if unicodedata.category(c) != 'Mn')


def select_mapping_type(mapping_types):
    mapping_order = ensembl_mapping_types
    selected_type = None
    mapping_types2 = [m[2] for m in mapping_types]
    for m in mapping_order:
        if m in mapping_types2:
            selected_type = m
            break
    assert selected_type is not None
    rv_eids = []
    for eid, canonical_name, mapping_type in mapping_types:
        if mapping_type == selected_type:
            rv_eids.append(eid)
    return selected_type, rv_eids


def contains_sublist(lst, sublst):
    if len(lst) < len(sublst):
        return False, None
    n = len(sublst)
    for i in range(len(lst) - n + 1):
        if sublst == lst[i:i + n]:
            return True, (i, i + n)
    return False, None


def eids_with_mapping_type(matches, mapping_type):
    rv_set = set()
    for (eid, canonical_name, m_mapping_type) in matches:
        if m_mapping_type == mapping_type:
            rv_set.add(eid)
    return list(rv_set)


def get_superword(start, end, indices):
    for other_start, other_end in indices:
        if start > other_start and end <= other_end:
            return other_start, other_end
        if start >= other_start and end < other_end:
            return other_start, other_end
    return None

alnum = re.compile('[a-zA-Z0-9]')
al = re.compile('[a-zA-Z]')


class GeneExtractor:
    def __init__(self, pubtator):
        self.gene_name_to_genes, self.gene_name_to_genes_permuted = gene_symbol_to_ensembl_id_map()
        self.pubtator = pubtator

    def extract(self, words, pmid):
        pubtator_mentions = self.pubtator.get_pmid2gene(pmid)
        all_words = words.split()
        rv = defaultdict(lambda: set())
        rv_mapping_types = {}
        rv_eids = defaultdict(lambda: set())
        for n in reversed(range(1, min(len(all_words), max_len))):
            for i in range(len(all_words) - n + 1):
                if not alnum.search(all_words[i]):
                    continue
                if not alnum.search(all_words[i + n - 1]):
                    continue
                words = [w for w in all_words[i:i + n]]
                orig_words = ' '.join(words)
                word_itself = normalize_gene_name(words, '')
                word_normalized = normalize_gene_name(words, 'PROTEIN')
                word_normalized_nospace = word_normalized.replace(' ', '')
                if not al.search(word_itself):
                    continue
                word = None
                if self.gene_name_to_genes_contains(word_itself, pubtator_mentions):
                    word = word_itself
                elif self.gene_name_to_genes_contains(word_normalized, pubtator_mentions):
                    word = word_normalized
                elif self.gene_name_to_genes_contains(word_normalized_nospace, pubtator_mentions):
                    word = word_normalized_nospace
                if word is not None:
                    matches = self.get_gene_name_to_genes_match(word, pubtator_mentions)
                    mapping_types = set()
                    for (eid, canonical_name, mapping_type) in matches:
                        mapping_types.add((eid, canonical_name, mapping_type))
                    mapping_type, eids = select_mapping_type(mapping_types)
                    if len(word) >= min_word_len[mapping_type]:
                        rv_mapping_types[word] = mapping_type
                        rv_eids[word].update(eids)
                        rv[word].add((i, i + n, orig_words))

        all_indices = set()
        for word in rv:
            all_indices.update(set((x[0], x[1]) for x in rv[word]))

        for word in rv:
            new_indices = set()
            new_orig_words = set()
            for start, end, orig_words in rv[word]:
                superword = get_superword(start, end, all_indices)
                if superword is None:
                    new_indices.add((start, end))
                else:
                    pass
                new_orig_words.add(orig_words)
            rv[word] = new_indices, new_orig_words

        result = defaultdict(set)        

        for word in rv:
            wordidxs, orig_words = rv[word]
            wordidxs_sort = sorted(list(wordidxs))
            wordidxs_start = [w[0] for w in wordidxs_sort]
            wordidxs_end = [w[1] for w in wordidxs_sort]
            mapping_type = rv_mapping_types[word]
            eids = rv_eids[word]
            surrounding_words_left = ['|^|'.join(all_words[start - 5:start]) for start in wordidxs_start]
            surrounding_words_right = ['|^|'.join(all_words[end:end + 5]) for end in wordidxs_end]
            gm = tuples.GeneMention(gene_name=word,
                                    wordidxs_start=tuple(sorted(list(wordidxs_start))),
                                    wordidxs_end=tuple(sorted(list(wordidxs_end))),
                                    num_mentions=len(rv[word][0]),
                                    mapping_type=mapping_type,
                                    surrounding_words_left=tuple(surrounding_words_left),
                                    surrounding_words_right=tuple(surrounding_words_right),
                                    orig_words='|^|'.join(list(orig_words)))
            for eid in eids:
                result[eid].add(gm)

        return result

    def gene_name_to_genes_contains(self, word, pubtator_mentions):
        if word in self.gene_name_to_genes:
            return True
        permuted = frozenset(word.split())
        if len(permuted) >= 3:
            if permuted in self.gene_name_to_genes_permuted:
                return True
        if word in pubtator_mentions:
            return True
        return False

    def get_gene_name_to_genes_match(self, word, pubtator_mentions):
        if word in pubtator_mentions:
            return pubtator_mentions[word]
        if word in self.gene_name_to_genes:
            return self.gene_name_to_genes[word]
        permuted = frozenset(word.split())
        if permuted in self.gene_name_to_genes_permuted:
            return self.gene_name_to_genes_permuted[permuted]
        return None
