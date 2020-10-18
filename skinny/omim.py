import urllib.request

import json
import re

import shutil
import io
import csv

import multiprocessing

from nltk.tokenize import sent_tokenize

from collections import defaultdict

import pkg_resources


def get_mim_names_to_ensemble():
    with urllib.request.urlopen('https://omim.org/static/omim/data/mim2gene.txt') as file:
        with io.TextIOWrapper(file) as text_file:
            for _ in range(4): text_file.readline() # Skip first 3 lines

            mim_to_ensemble = {}

            for row in csv.DictReader(text_file, delimiter='\t'):
                if row['Ensembl Gene ID (Ensembl)'] == '':
                    continue

                mim_to_ensemble[int(row['# MIM Number'])] = row['Ensembl Gene ID (Ensembl)']

            return mim_to_ensemble


def get_omin_file(mim_number, api_key, cache_directory=None):
    path = 'https://api.omim.org/api/entry?mimNumber={}&include=all&apiKey={}&format=json'.format(mim_number, api_key)

    last_exception = None
    # 5 attempts
    for i in range(5):
        try:
            with urllib.request.urlopen(path.format(api_key)) as file:
                with io.TextIOWrapper(file) as text_file:
                    omim_entry = json.load(text_file)
                    return omim_entry["omim"]["entryList"][0]["entry"]
        except Exception as e:
            last_exception = e
    else:
        if last_exception is not None:
            raise last_exception


variant_reference_regex = re.compile(r'{([0-9]+):')


def get_references_from_variant_text(variant_text):
    variant_text = re.sub('((previously|originally) (described|reported|studied) by [^{]*{[^}]*}[^,\.]*[\.,])',
                          '', variant_text)
    locations = variant_reference_regex.finditer(variant_text)
    return {int(location.group(1)) for location in locations}


def get_disease_database(omim_map, hgmd_location, mim_to_ensemble):
    pmid_to_info = defaultdict(lambda: {
        'variant_types': set(),
        'inheritance_modes': set(),
        'topical_genes': set(),
    })

    for name, entry in omim_map.items():
        if 'allelicVariantList' not in entry:
            continue

        if 'referenceList' not in entry:
            continue

        pmids_for_reference = {}

        for reference in entry['referenceList']:
            reference_entry = reference['reference']
            if 'pubmedID' in reference_entry:
                pmids_for_reference[reference_entry['referenceNumber']] = reference_entry['pubmedID']

        for variant in entry['allelicVariantList']:
            variant_entry = variant['allelicVariant']
            if 'text' in variant_entry:
                paragraphs = variant_entry['text'].split('\n\n')

                topical_gene = mim_to_ensemble[name]

                if 'mutations' in variant_entry:
                    variant_types = set(get_variant_type(variant_entry['mutations']))
                else:
                    variant_types = set()

                for para in paragraphs:
                    references_for_para = get_references_from_variant_text(para)

                    inheritance_modes = get_inheritance_modes(para)

                    for reference in references_for_para:
                        if reference in pmids_for_reference:
                            pmid = pmids_for_reference[reference]

                            pmid_to_info[pmid]['variant_types'].update(variant_types)
                            pmid_to_info[pmid]['inheritance_modes'].update(inheritance_modes)
                            pmid_to_info[pmid]['topical_genes'].add(topical_gene)

    with open(hgmd_location, 'r') as hgmd:
        reader = csv.DictReader(hgmd, delimiter='\t')

        for row in reader:
            if row['tag'] != 'DM':
                continue

            try:
                pmid = int(row['pmid'])
            except ValueError:
                continue

            try:
                mim_id = int(row['omimid'])
            except ValueError:
                continue

            if mim_id in mim_to_ensemble:
                topical_gene = mim_to_ensemble[mim_id]
                pmid_to_info[pmid]['topical_genes'].add(topical_gene)
            else:
                pmid_to_info[pmid] # Force the defaultdict to have this

    # Convert to lists for JSON serialization
    for pmid, data in pmid_to_info.items():
        for key in data:
            data[key] = list(data[key])
   
    return pmid_to_info


def get_inheritance_modes(paragraph_text):
    dominant = False
    recessive = False
    sentences = sent_tokenize(paragraph_text)
    for sent in sentences:
        if 'parent' in sent or 'brother' in sent or 'sister' in sent or 'sibling' in sent or 'family' in sent or 'mother' in sent or 'father' in sent:
            continue
        sent = sent.split()
        for i, word in enumerate(sent):
            word = word.lower()
            if word == 'homozygous' or word == 'homozygote' or word == 'homozygosity' or word == 'recessive':
                recessive = True
            if len(sent) > i+1 and word == 'compound' and sent[i+1].startswith('het'):
                recessive = True
            if word == 'heterozygous' or word == 'heterozygote' or word == 'heterozygosity':
                if i > 0:
                    if sent[i-1] != 'compound':
                        dominant = True
                else:
                    dominant = True
            if word == 'dominant':
                dominant = True

    result = set()

    if dominant:
        result.add('dominant')

    if recessive:
        result.add('recessive')

    return result


lof = "destructive"
gof = "not_destructive"

aa = "(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR)"

missense_pattern = re.compile('%s[0-9]+%s' % (aa, aa))

stopgain_pattern = re.compile('%s[0-9]+TER' % aa)
stoploss_pattern = re.compile('TER[0-9]+%s' % aa)
splicing_pattern = re.compile('IVS[0-9]+')

del_pattern = re.compile('([0-9]+)-BP DEL')
dup_pattern = re.compile('([0-9]+)-BP DUP')
ins_pattern = re.compile('([0-9]+)-BP INS')

gof_patterns = [missense_pattern]
lof_patterns = [stopgain_pattern, stoploss_pattern, splicing_pattern]
lof_or_gof_patterns = [del_pattern, dup_pattern, ins_pattern]


def get_variant_type(mutation):
    def lof_or_gof(m):
        num_deleted = int(m.group(1))
        if num_deleted % 3 == 0:
            if num_deleted <= 9:
                return gof
            else:
                return lof
        else:
            return lof

    for p in gof_patterns:
        m = p.search(mutation)
        if m:
            yield gof
    for p in lof_patterns:
        m = p.search(mutation)
        if m:
            yield lof
    for p in lof_or_gof_patterns:
        m = p.search(mutation)
        if m:
            yield lof_or_gof(m)
