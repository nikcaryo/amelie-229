
import re
import nltk
nltk.data.path.append('/cluster/u/jbirgmei/nltk_data')
nltk.download('wordnet')
from nltk.stem.wordnet import WordNetLemmatizer
nltk.download('stopwords')
from nltk.corpus import stopwords

from .dataloaders.hpo import HPO
from . import tuples


# some of this junk is never used:
max_len = 8
min_word_len = 3
split_list = [',', ';', '.']
split_max_stops = 2
permuted = True
omitted_interior = True
VALS = [('neg', False), ('pos', True)]

lmtzr = WordNetLemmatizer()


def lemmatize(w):
    if w.isalpha():
        return lmtzr.lemmatize(w)
    else:
        return w


def normalize_phrase(p):
    """Lowercases, removes stop words, and lemmatizes inputted multi-word phrase"""
    out = []

    # split into contiguous alphanumeric segments, lower-case, filter stopwords, lemmatize
    ws = [re.sub(r'[^a-z0-9]', '', w) for w in p.lower().split()]
    ws = [w for w in ws if w not in stopwords.words('english')]
    ws = [lemmatize(w) for w in ws]
    out.append(' '.join(ws))

    # if there's a comma, try permuting the order (some of these cases omitted from HPO!)
    if ',' in p:
        cs = re.split(r'\s*,\s*', p.strip())
        out += normalize_phrase(' '.join(cs[::-1]))
    return out


forbidden_forms = ['cell']


def enrich_pheno(form):
    ret = []
    hpoid, phrase, entry_type, canon_name = [x.strip() for x in form]
    if 'decreased' in phrase:
        new_pheno = phrase.replace('decreased', 'reduced')
        ret.append((hpoid, new_pheno, 'MORPHED', canon_name))

    new_ret = []
    for row in ret:
        hpoid, pheno, entry_type, canon_name = [x.strip() for x in row]
        words = pheno.split()
        for word in words:
            # just assuming that only one slash occurs per line
            if '/' in word:
                nword = []
                nword.append(word.split('/')[0])
                nword.append(word.split('/')[1])
                new_pheno = pheno.replace(word, nword[0])
                new_ret.append((hpoid, new_pheno, 'SLASHED', canon_name))
                new_pheno = pheno.replace(word, nword[1])
                new_ret.append((hpoid, new_pheno, 'SLASHED', canon_name))
    ret = ret + new_ret
    return ret


lemmatizer = WordNetLemmatizer()


def process_pheno_terms(hpo, xref_synonyms=False):
    seen = {}
    rv = []
    for term in hpo.terms:
        hpo_id = term
        canon_name = hpo.terms[term]['name'][0].strip()
        exact = [hpo.terms[term]['name'][0].strip().lower()]
        for synonym in hpo.get_synonyms(term):
            exact.append(synonym.replace('"', ' ').strip().lower())
        if xref_synonyms:
            for synonym in hpo.get_xref_synonyms(term):
                exact.append(synonym.replace('"', ' ').strip().lower())
        forms = [(hpo_id, p, "EXACT", canon_name) for p in exact]
        for p in exact:
            forms += [(hpo_id, np, "LEMMA", canon_name) for np in normalize_phrase(p) if len(np.strip()) > 0]
            enriched = enrich_pheno((hpo_id, p, "EXACT", canon_name))
            forms.extend(enriched)
            for e in enriched:
                forms += [(hpo_id, np, "%s_LEMMA" % e[2], canon_name)
                          for np in normalize_phrase(e[1]) if len(np.strip()) > 0]
        for f in forms:
            k = f[0] + f[1]
            if k not in seen:
                seen[k] = 1
                if f[1] not in forbidden_forms:
                    rv.append(f)
    return rv


def load_pheno_terms(hpo, xref_synonyms=False):
    phenos = {}
    pheno_sets = {}
    rows = process_pheno_terms(hpo, xref_synonyms=xref_synonyms)
    for row in rows:
        # print(row)
        hpoid, phrase, entry_type, canon_name = [x.strip() for x in row]
        if hpoid in hpo.terms:
            if phrase in phenos:
                phenos[phrase].append((hpoid, entry_type, canon_name))
            else:
                phenos[phrase] = [(hpoid, entry_type, canon_name)]
            phrase_bow = frozenset(phrase.split())
            if phrase_bow in pheno_sets:
                pheno_sets[phrase_bow].append((hpoid, entry_type, canon_name))
            else:
                pheno_sets[phrase_bow] = [(hpoid, entry_type, canon_name)]
    return phenos, pheno_sets


def keep_word(w):
    return w.lower() not in stopwords.words("english") and len(w) >= 3


reW = re.compile('\W+')


class PhenoExtractor:

    def __init__(self, xref_synonyms=False):
        self.hpo = HPO(xref_synonyms=xref_synonyms)
        phenos, pheno_sets = load_pheno_terms(self.hpo, xref_synonyms=xref_synonyms)
        self.terms = phenos
        self.term_sets = pheno_sets

    def create_supervised_mention(self, all_words, idxs, entity, canon_name):
        words = [all_words[i] for i in idxs]
        pm = tuples.PhenoMention(wordidxs=idxs, words=words, hpo_id=entity, canon_name=canon_name)
        return pm

    def extract(self, words):
        all_words = words.split()
        all_lemmas = [lmtzr.lemmatize(w.lower()) for w in all_words]
        split_indices = set()
        split_indices.update([i for i, w in enumerate(all_words) if w in split_list])

        seq = []
        for i,w in enumerate(all_words):
            if not keep_word(w):
                seq.append(i)
            else:
                if len(seq) > split_max_stops:
                    split_indices.update(seq)
                seq = []

        for n in reversed(range(1, min(len(all_words), max_len))):
            for i in range(len(all_words)-n+1):
                wordidxs = [j for j in range(i,i+n)]
                words = [reW.sub(' ', w.lower()) for w in all_words[i:i+n]]
                lemmas = [reW.sub(' ', w) for w in all_lemmas[i:i+n]]

                if not split_indices.isdisjoint(wordidxs):
                    continue
                if not all(map(keep_word, [words[0], lemmas[0], words[-1], lemmas[-1]])):
                    continue
                ws, lws = zip(*[(words[k], lemmas[k]) for k in range(n) if keep_word(words[k]) and keep_word(lemmas[k])])
                p, lp = map(' '.join, [ws, lws])

                if p in self.terms or lp in self.terms:
                    entities = self.terms[p] if p in self.terms else self.terms[lp]
                    for (entity, entry_type, canon_name) in entities:
                        m = self.create_supervised_mention(all_words, wordidxs, 
                                                           entity, canon_name)
                        if m is not None:
                            yield m
                    split_indices.update(wordidxs)
                    continue

                ps, lps = map(frozenset, [ws, lws])
                if (len(ps) == len(ws) and ps in self.term_sets) or (len(lps) == len(lws) and lps in self.term_sets):
                    entities = self.term_sets[ps] if ps in self.term_sets else self.term_sets[lps]
                    for (entity, entry_type, canon_name) in entities:
                        m = self.create_supervised_mention(all_words, wordidxs, 
                                                           entity, canon_name)
                        if m is not None:
                            yield m
                    continue

                if len(ws) > 2:
                    for omit in range(1, len(ws)-1):
                        p, lp = [' '.join([w for i,w in enumerate(x) if i != omit]) for x in [ws, lws]]
                        if p in self.terms or lp in self.terms:
                            entities = self.terms[p] if p in self.terms else self.terms[lp]
                            for (entity, entry_type, canon_name) in entities:
                                m = self.create_supervised_mention(all_words, wordidxs, 
                                                                   entity, canon_name)
                                if m is not None:
                                    yield m


def num_pheno_terms():
    hpo = HPO(xref_synonyms=False)
    rv = process_pheno_terms(hpo, xref_synonyms=False)
    hpo_ids = set()
    phrases = set()
    for hpoid, phrase, entry_type, canon_name in rv:
        hpo_ids.add(hpoid)
        phrases.add(phrase)
    print("Non XREF: %d HPO IDs, %d phrases" % (len(hpo_ids), len(phrases)))

    hpo = HPO(xref_synonyms=True)
    rv = process_pheno_terms(hpo, xref_synonyms=True)
    hpo_ids = set()
    phrases = set()
    for hpoid, phrase, entry_type, canon_name in rv:
        hpo_ids.add(hpoid)
        phrases.add(phrase)
    print("With XREF: %d HPO IDs, %d phrases" % (len(hpo_ids), len(phrases)))

