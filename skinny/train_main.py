import argparse
import json
import multiprocessing
import os
import pickle
import random
import sys
from collections import defaultdict

import constdb

from amelie.dataloaders import clinvar
from . import pheno_extractor, gene_extractor, process_main
from .dataloaders import pubtator
from . import text_classification, topical_gene, variant_classification
from .tuples import VCFCoords
from amelie.amelie.fake_patients_main import load_name_to_ensembl


def load_article(processed_dir, pmid):
    path = processed_dir + '/' + str(pmid) + '.pkl'

    if not os.path.isfile(path):
        return None
    else:
        with open(path, 'rb') as file:
            try:
                loaded = pickle.load(file)
                return loaded
            except EOFError:
                print("Cannot load pmid %s for some reason ..." % pmid, flush=True)
                return None


class ConvertToTextFunction:
    def __init__(self, processed_dir, replace_phenos_with_nothing):
        self.processed_dir = processed_dir
        self.replace_phenos_with_nothing = replace_phenos_with_nothing

    def __call__(self, pmid):
        processed_article = load_article(self.processed_dir, pmid)
        if processed_article is None:
            return None

        article = text_classification.convert_to_text(processed_article, use_main_text=True,
                                                      replace_phenos_with_nothing=self.replace_phenos_with_nothing)

        return pmid, article


def convert_all_to_text(processed_dir, pmids, replace_phenos_with_nothing):
    with multiprocessing.Pool(100) as pool:
        print('Pool created', flush=True)
        for i, item in enumerate(pool.imap_unordered(
                ConvertToTextFunction(processed_dir,
                                      replace_phenos_with_nothing=replace_phenos_with_nothing),
                pmids,
                chunksize=100)):
            if i % 10000 == 0:
                print('Processed ', i, ' out of ', len(pmids), flush=True)

            if item is not None:
                yield item


class ConvertHeaderAbstractToTextFunction:
    def __init__(self, out_dir, replace_phenos_with_nothing, process_dir, ta_genephenos_dir):
        self.out_dir = out_dir
        self.replace_phenos_with_nothing = replace_phenos_with_nothing
        self.process_dir = process_dir
        self.ta_genephenos_dir = ta_genephenos_dir

    def __call__(self, pubmed_entries):
        # keys = {pmid for pmid,_ in pubmed_entries}
        print("ConvertHeaderAbstractToTextFunction", flush=True)
        with pubtator.Pubtator(self.out_dir + '/pubtator.db') as pub:
            print("Loaded PubTator", flush=True)
            extractors = {
                'p_extractor': pheno_extractor.PhenoExtractor(),
                'g_extractor': gene_extractor.GeneExtractor(pub)
            }
            print("Loaded gene and pheno extractor", flush=True)

            results = []

            for i, (pmid, pubmed_info) in enumerate(pubmed_entries):
                print(i, flush=True)
                processed_file = os.path.join(self.process_dir, '%s.pkl' % pmid)
                ta_genephenos_file = os.path.join(self.ta_genephenos_dir, '%s.pkl' % pmid)
                if os.path.isfile(processed_file):
                    with open(processed_file, 'rb') as f:
                        text_results = pickle.load(f)
                elif os.path.isfile(ta_genephenos_file):
                    with open(ta_genephenos_file, 'rb') as f:
                        text_results = pickle.load(f)
                else:
                    text_results = {
                        'abstract': process_main.process_part_of_paper('abstract', pubmed_info['abstract'], extractors, pmid),
                        'title': process_main.process_part_of_paper('title', pubmed_info['title'], extractors, pmid),
                    }
                    with open(ta_genephenos_file, 'wb') as f:
                        pickle.dump(text_results, f)
                article = text_classification.\
                    convert_to_text(text_results, use_main_text=False,
                                    replace_phenos_with_nothing=self.replace_phenos_with_nothing)
                results.append((pmid, article))

            return results


def convert_header_abstract_to_text(pubmed_entries, out_dir, replace_phenos_with_nothing,
                                    process_dir, ta_genephenos_dir):
    print("Converting header abstract to text", flush=True)

    num_processes = 20

    num_items_per_batch = (len(pubmed_entries) + num_processes - 1) // num_processes

    batches = [pubmed_entries[i:i+num_items_per_batch] for i in range(0, len(pubmed_entries), num_items_per_batch)]

    print(num_processes, flush=True)
    print(len(batches), flush=True)

    with multiprocessing.Pool(num_processes) as pool:
        print('Pool created', flush=True)
        for i, items in enumerate(pool.imap_unordered(
                ConvertHeaderAbstractToTextFunction(out_dir,
                                                    replace_phenos_with_nothing=replace_phenos_with_nothing,
                                                    process_dir=process_dir,
                                                    ta_genephenos_dir=ta_genephenos_dir),
                batches,
                chunksize=1)):
            yield from items


class ConvertToArticleFunction:
    def __init__(self, processed_dir):
        self.processed_dir = processed_dir

    def __call__(self, pmid):
        processed_article = load_article(self.processed_dir, pmid)
        if processed_article is None:
            return None

        return pmid, processed_article


def convert_all_to_articles(processed_dir, pmids):
    with multiprocessing.Pool(100) as pool:
        print('Pool created', flush=True)
        for i, item in enumerate(pool.imap_unordered(
                ConvertToArticleFunction(processed_dir),
                pmids,
                chunksize=100)):
            if i % 10000 == 0:
                print('Processed ', i, ' out of ', len(pmids), flush=True)

            if item is not None:
                yield item


class ConvertToVariantsFunction:
    def __init__(self, processed_dir, clinvar_db, topical_gene_classifier):
        self.processed_dir = processed_dir
        self.clinvar_db = clinvar_db
        self.topical_gene_classifier = topical_gene_classifier

    def __call__(self, pmid):
        def modification_to_vcf_coords(m):
            assert m.vcf_pos is None or type(m.vcf_pos) == int
            return VCFCoords(chrom=m.chrom.lstrip('chr'),
                             pos=m.vcf_pos,
                             ref=m.vcf_ref,
                             alt=m.vcf_alt)

        processed_article = load_article(self.processed_dir, pmid)
        if processed_article is None:
            return None

        expected_changes = self.clinvar_db.get_dna_changes(pmid)
        variants = processed_article['text']['variant_mentions']
        if len(variants) == 0:
            return None
        variants_by_text = defaultdict(set)
        for variant_info in variants:
            variants_by_text[variant_info.variant_text].add(variant_info)
        extractor = variant_classification.DocumentVariantFeatureExtractor(processed_article,
                                                                           self.topical_gene_classifier)
        results = []

        for original_variant_string, variants_info in variants_by_text.items():
            has_one_positive = False
            for variant_info in variants_info:
                is_valid_variant = any(modification_to_vcf_coords(modification) in expected_changes
                                       for modification in variant_info.modifications)
                if is_valid_variant:
                    has_one_positive = True

            if has_one_positive:
                for variant_info in variants_info:
                    is_valid_variant = any(modification_to_vcf_coords(modification) in expected_changes
                                           for modification in variant_info.modifications)
                    overlapping_modifications = []
                    for modification in variant_info.modifications:
                        if modification_to_vcf_coords(modification) in expected_changes:
                            overlapping_modifications.append(modification)
                    variant_rows = extractor.featurize(variant_info)
                    for row in variant_rows:
                        results.append((is_valid_variant, row))

        return results


def convert_all_to_variants(processed_dir, pmids, clinvar_db, topical_gene_classifier):
    # pmids = pmids[:500]  # for testing
    with multiprocessing.Pool(40) as pool:
        print('Pool created', flush=True)
        for i, item in enumerate(pool.imap_unordered(
                ConvertToVariantsFunction(processed_dir, clinvar_db, topical_gene_classifier),
                pmids,
                chunksize=100
            )):
            if i % 1000 == 0:
                print('Processed ', i, ' out of ', len(pmids), flush=True)

            if item is not None:
                for sub_item in item:
                    yield sub_item


def train_relevance_classifier(out_dir, processed_dir, allowed_pmids, forbidden_pmids,
                               save, cross_val, l1, replace_phenos_with_nothing):
    with open(out_dir + '/dataset_meta.json') as file:
        dataset_info = json.load(file)

    positives = set(dataset_info['positive_pmids']) - forbidden_pmids
    negatives = set(dataset_info['negative_pmids']) - forbidden_pmids
    print('Total positives: ', len(positives), flush=True)
    print('Total negatives: ', len(negatives), flush=True)

    print('Starting to pull article information', flush=True)

    articles = []
    labels = []
    if allowed_pmids is not None:
        total_to_process = list(set(positives | negatives).intersection(allowed_pmids))
        print("Skipping %d PMIDs because not in allowed PMIDs" % (len(list(set(positives | negatives)))
                                                                  - len(total_to_process)), flush=True)
    else:
        total_to_process = list(positives | negatives)

    for pmid, article in convert_all_to_text(processed_dir, total_to_process,
                                             replace_phenos_with_nothing=replace_phenos_with_nothing):
        label = pmid in positives
        articles.append(article)
        labels.append(label)

    print('Have ', len(articles), flush=True)
    print('Num positives ', sum(1 if a else 0 for a in labels), flush=True)
    print('Num negatives ', sum(1 if not a else 0 for a in labels), flush=True)
    print('Done converting', flush=True)

    classifier = text_classification.create_model(articles, labels, cross_val=cross_val, l1=l1)
    if save:
        print("SAVING PDF RELEVANCE CLASSIFIER", flush=True)
        with open(out_dir + '/pdf_relevance.pkl', 'wb') as out_file:
            pickle.dump(classifier, out_file)
    else:
        print("NOT SAVING PDF RELEVANCE CLASSIFIER", flush=True)


def train_title_abstract_relevance_classifier(out_dir, process_dir, ta_genephenos_dir,
                                              allowed_pmids=None, forbidden_pmids=None,
                                              save=True, cross_val=False, l1=False,
                                              replace_phenos_with_nothing=False,
                                              use_processed=True):
    with open(out_dir + '/dataset_meta.json') as file:
        dataset_info = json.load(file)

    positives = set(dataset_info['positive_pmids']) - forbidden_pmids

    with constdb.read(out_dir + '/pubmed.db') as pubmed_dump:
        pubmed_pmids = set(pubmed_dump.keys())

    negatives = pubmed_pmids - positives - forbidden_pmids

    if allowed_pmids is not None:
        orig_len_positives = len(positives)
        orig_len_negatives = len(negatives)
        positives = positives.intersection(allowed_pmids)
        negatives = negatives.intersection(allowed_pmids)
        print("Skipping %d positives b/c of allowed PMIDs" % (orig_len_positives - len(positives)), flush=True)
        print("Skipping %d negatives b/c of allowed PMIDs" % (orig_len_negatives - len(negatives)), flush=True)
    print('Total positives: ', len(positives), flush=True)
    print('Total negatives: ', len(negatives), flush=True)

    negative_ratio = 1.2
    negatives = set(random.sample(list(negatives), int(len(positives) * negative_ratio)))
    print('Sampled negatives: ', len(negatives), flush=True)

    print('Starting to pull article information', flush=True)

    articles = []
    labels = []

    total_to_process = list(positives | negatives)

    entries = []

    # if cross_val:
    #     random.shuffle(total_to_process)
    #     print("Cross-validating TA classifier; downsampling to 5,000")
    #     total_to_process = total_to_process[:5000]

    with constdb.read(out_dir + '/pubmed.db', keys_to_read=set(total_to_process)) as pubmed_dump:
        print('Pubmed dump opened', flush=True)
        for i, pmid in enumerate(total_to_process):
            if (i + 1) % 1000 == 0:
                print('Done with', i, 'out of', len(total_to_process), flush=True)
            entries.append((pmid, json.loads(pubmed_dump.get(pmid))))

    print('Starting to convert data', flush=True)

    for i, (pmid, article) in enumerate(
            convert_header_abstract_to_text(entries, out_dir,
                                            replace_phenos_with_nothing=replace_phenos_with_nothing,
                                            process_dir=process_dir,
                                            ta_genephenos_dir=ta_genephenos_dir)):
        print(i, flush=True)
        label = pmid in positives
        articles.append(article)
        labels.append(label)

    print('Have ', len(articles), flush=True)
    print('Num positives', sum(1 if a else 0 for a in labels), flush=True)
    print('Num negatives', sum(1 if not a else 0 for a in labels), flush=True)
    print('Done converting', flush=True)

    classifier = text_classification.create_model(articles, labels, cross_val=cross_val, l1=l1)
    if save:
        print("SAVING PUBMED ONLY RELEVANCE CLASSIFIER", flush=True)
        with open(out_dir + '/pubmed_only_relevance.pkl', 'wb') as out_file:
            pickle.dump(classifier, out_file)
    else:
        print("NOT SAVING PUBMED ONLY RELEVANCE CLASSIFIER", flush=True)


def process_train_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('out_dir', type=str,
                        default='pubmed_relevance_classifier.pkl',
                        help='The name of the database to output')

    parser.add_argument('process_dir', type=str,
                        default='dataset_meta.json',
                        help='The name of the database to output')

    parser.add_argument('ta_genephenos_dir', type=str,
                        default='dataset_meta.json',
                        help='The name of the database to output')

    parser.add_argument('--article_max_year', type=int, default=None)
    parser.add_argument('--article_to_year_map', type=str, default=None)
    parser.add_argument('--forbidden_genesymbols_file', type=str, default=None)
    parser.add_argument('--title_abstract_only', action='store_true', default=False)
    parser.add_argument('--do_not_save', action='store_true', default=False)
    parser.add_argument('--cross_validate', action='store_true', default=False)
    parser.add_argument('--l1_regularization', action='store_true', default=False)
    parser.add_argument('--replace_phenos_with_nothing', action='store_true', default=False)

    args = parser.parse_args()

    if not os.path.isdir(args.out_dir):
        print("Out dir %s does not exist" % args.out_dir, flush=True)
        sys.exit(1)

    if not os.path.isdir(args.process_dir):
        print("Processed dir %s does not exist" % args.process_dir, flush=True)
        sys.exit(1)

    return args


def load_allowed_pmids(max_year, pmid_to_year_map):
    if max_year is None:
        return None
    rv = set()
    with open(pmid_to_year_map) as f:
        for line in f:
            line = line.rstrip('\n')
            pmid, year = line.split('\t')
            year = int(year)
            if year <= max_year:
                rv.add(pmid)
                rv.add(int(pmid))
    return rv


def load_forbidden_pmids(forbidden_genesymbols_filename, out_dir):
    if forbidden_genesymbols_filename is None:
        return set()
    with open(forbidden_genesymbols_filename) as f:
        forbidden_gene_symbols = set([x.strip() for x in f.readlines() if len(x.strip()) > 0])
    name_to_ensembl = load_name_to_ensembl()
    forbidden_ensgenes = set()
    for gene_symbol in forbidden_gene_symbols:
        assert gene_symbol in name_to_ensembl, gene_symbol
        forbidden_ensgenes.add(name_to_ensembl[gene_symbol])
    with open(out_dir + '/omim.json') as file:
        omim = json.load(file)
    rv = set()
    for pmid in omim:
        topical_genes = set(omim[pmid]["topical_genes"])
        if len(topical_genes.intersection(forbidden_ensgenes)) > 0:
            rv.add(str(pmid))
            rv.add(int(pmid))
    print("Have %d forbidden pmids (%s, ...)" % (len(rv), ', '.join(list(str(x) for x in rv)[:10])), flush=True)
    return rv


def train_classifiers_program():
    args = process_train_args()

    allowed_pmids = load_allowed_pmids(max_year=args.article_max_year, pmid_to_year_map=args.article_to_year_map)
    print("Loading forbidden PMIDs", flush=True)
    forbidden_pmids = load_forbidden_pmids(args.forbidden_genesymbols_file, out_dir=args.out_dir)
    print("Training title/abstract relevance classifier", flush=True)
    train_title_abstract_relevance_classifier(args.out_dir,
                                              args.process_dir,
                                              args.ta_genephenos_dir,
                                              allowed_pmids=allowed_pmids,
                                              forbidden_pmids=forbidden_pmids,
                                              save=not args.do_not_save,
                                              cross_val=args.cross_validate,
                                              l1=args.l1_regularization,
                                              replace_phenos_with_nothing=args.replace_phenos_with_nothing)
    print("Training full text relevance classifier", flush=True)
    train_relevance_classifier(args.out_dir, args.process_dir, allowed_pmids=allowed_pmids,
                               forbidden_pmids=forbidden_pmids,
                               save=not args.do_not_save,
                               cross_val=args.cross_validate,
                               l1=args.l1_regularization,
                               replace_phenos_with_nothing=args.replace_phenos_with_nothing)
    print("Training inheritance mode classifier", flush=True)
    train_text_field_classifier(args.out_dir, args.process_dir, 'inheritance_modes',
                                allowed_pmids=allowed_pmids,
                                forbidden_pmids=forbidden_pmids,
                                save=not args.do_not_save,
                                cross_val=args.cross_validate,
                                l1=args.l1_regularization,
                                replace_phenos_with_nothing=args.replace_phenos_with_nothing)
    print("Training variant type classifier", flush=True)
    train_text_field_classifier(args.out_dir, args.process_dir, 'variant_types',
                                allowed_pmids=allowed_pmids,
                                forbidden_pmids=forbidden_pmids,
                                save=not args.do_not_save,
                                cross_val=args.cross_validate,
                                l1=args.l1_regularization,
                                replace_phenos_with_nothing=args.replace_phenos_with_nothing)
    print("Training topical gene classifier", flush=True)
    train_topical_gene_classifier(args.out_dir, args.process_dir,
                                  allowed_pmids=allowed_pmids,
                                  forbidden_pmids=forbidden_pmids,
                                  save=not args.do_not_save,
                                  cross_val=args.cross_validate,
                                  l1=args.l1_regularization,
                                  replace_phenos_with_nothing=args.replace_phenos_with_nothing)
    print("Training variant classifier", flush=True)
    train_variant_classifier(args.out_dir, args.process_dir,
                             allowed_pmids=allowed_pmids,
                             forbidden_pmids=forbidden_pmids,
                             save=not args.do_not_save,
                             cross_val=args.cross_validate,
                             replace_phenos_with_nothing=args.replace_phenos_with_nothing)


def train_title_abstract():
    args = process_train_args()
    print("Training title/abstract relevance classifier", flush=True)
    allowed_pmids = load_allowed_pmids(max_year=args.article_max_year, pmid_to_year_map=args.article_to_year_map)
    print("Loading forbidden PMIDs", flush=True)
    forbidden_pmids = load_forbidden_pmids(args.forbidden_genesymbols_file, out_dir=args.out_dir)
    train_title_abstract_relevance_classifier(args.out_dir,
                                              args.process_dir,
                                              args.ta_genephenos_dir,
                                              allowed_pmids=allowed_pmids,
                                              forbidden_pmids=forbidden_pmids,
                                              save=not args.do_not_save,
                                              cross_val=args.cross_validate,
                                              l1=args.l1_regularization,
                                              replace_phenos_with_nothing=args.replace_phenos_with_nothing)


def train_full_text_relevance():
    args = process_train_args()
    print("Training full text relevance classifier", flush=True)
    allowed_pmids = load_allowed_pmids(max_year=args.article_max_year, pmid_to_year_map=args.article_to_year_map)
    print("Loading forbidden PMIDs", flush=True)
    forbidden_pmids = load_forbidden_pmids(args.forbidden_genesymbols_file, out_dir=args.out_dir)
    train_relevance_classifier(args.out_dir, args.process_dir,
                               allowed_pmids=allowed_pmids,
                               forbidden_pmids=forbidden_pmids,
                               save=not args.do_not_save,
                               cross_val=args.cross_validate,
                               l1=args.l1_regularization,
                               replace_phenos_with_nothing=args.replace_phenos_with_nothing)


def train_inheritance_mode():
    args = process_train_args()
    print("Training inheritance mode classifier", flush=True)
    allowed_pmids = load_allowed_pmids(max_year=args.article_max_year, pmid_to_year_map=args.article_to_year_map)
    print("Loading forbidden PMIDs", flush=True)
    forbidden_pmids = load_forbidden_pmids(args.forbidden_genesymbols_file, out_dir=args.out_dir)
    train_text_field_classifier(args.out_dir, args.process_dir, 'inheritance_modes',
                                allowed_pmids=allowed_pmids,
                                forbidden_pmids=forbidden_pmids,
                                save=not args.do_not_save,
                                cross_val=args.cross_validate,
                                l1=args.l1_regularization,
                                replace_phenos_with_nothing=args.replace_phenos_with_nothing)


def train_variant_type():
    args = process_train_args()
    print("Training variant type classifier", flush=True)
    allowed_pmids = load_allowed_pmids(max_year=args.article_max_year, pmid_to_year_map=args.article_to_year_map)
    forbidden_pmids = load_forbidden_pmids(args.forbidden_genesymbols_file, out_dir=args.out_dir)
    train_text_field_classifier(args.out_dir, args.process_dir, 'variant_types',
                                allowed_pmids=allowed_pmids,
                                forbidden_pmids=forbidden_pmids,
                                save=not args.do_not_save,
                                cross_val=args.cross_validate,
                                l1=args.l1_regularization,
                                replace_phenos_with_nothing=args.replace_phenos_with_nothing)


def train_topical_gene():
    args = process_train_args()
    print("Training topical gene classifier", flush=True)
    allowed_pmids = load_allowed_pmids(max_year=args.article_max_year, pmid_to_year_map=args.article_to_year_map)
    print("Loading forbidden PMIDs", flush=True)
    forbidden_pmids = load_forbidden_pmids(args.forbidden_genesymbols_file, out_dir=args.out_dir)
    train_topical_gene_classifier(args.out_dir, args.process_dir,
                                  allowed_pmids=allowed_pmids,
                                  forbidden_pmids=forbidden_pmids,
                                  save=not args.do_not_save,
                                  cross_val=args.cross_validate,
                                  l1=args.l1_regularization,
                                  replace_phenos_with_nothing=args.replace_phenos_with_nothing)


def train_variants():
    args = process_train_args()
    print("Training variant classifier", flush=True)
    allowed_pmids = load_allowed_pmids(max_year=args.article_max_year, pmid_to_year_map=args.article_to_year_map)
    print("Loading forbidden PMIDs", flush=True)
    forbidden_pmids = load_forbidden_pmids(args.forbidden_genesymbols_file, out_dir=args.out_dir)
    train_variant_classifier(args.out_dir, args.process_dir,
                             allowed_pmids=allowed_pmids,
                             forbidden_pmids=forbidden_pmids,
                             save=not args.do_not_save,
                             cross_val=args.cross_validate,
                             replace_phenos_with_nothing=args.replace_phenos_with_nothing)


def shuffle_articles_labels(articles, labels):
    zipped_articles_labels = [x for x in zip(articles, labels)]
    random.shuffle(zipped_articles_labels)
    return [x for x in zip(*zipped_articles_labels)]


def train_text_field_classifier(out_dir, process_dir, field_name, allowed_pmids,
                                forbidden_pmids, save, cross_val, l1,
                                replace_phenos_with_nothing):
    with open(out_dir + '/dataset_meta.json') as file:
        dataset_info = json.load(file)

    with open(out_dir + '/omim.json') as file:
        omim = json.load(file)
    positives = set(dataset_info['positive_pmids']) - forbidden_pmids
    print('Total positives: ', len(positives), flush=True)
    good_positives = []

    for pmid in positives:
        omim_data = omim[str(pmid)]
        if field_name not in omim_data:
            continue
        field = omim_data[field_name]
        if len(field) == 0 or len(field) == 2:
            # Skip bad ones
            continue
        good_positives.append(pmid)

    articles = []
    labels = []

    for pmid, article in convert_all_to_text(process_dir, good_positives,
                                             replace_phenos_with_nothing=replace_phenos_with_nothing):
        if allowed_pmids is not None and pmid not in allowed_pmids:
            print("Skipping PMID %s because not in allowed PMIDs" % pmid, flush=True)
            continue
        omim_data = omim[str(pmid)]
        field = omim_data[field_name]
        articles.append(article)
        labels.append(field[0])

    for pmid, article in convert_all_to_text(process_dir, dataset_info['gwas_pmids'],
                                             replace_phenos_with_nothing=replace_phenos_with_nothing):
        if allowed_pmids is not None and pmid not in allowed_pmids:
            print("Skipping PMID %s because not in allowed PMIDs" % pmid, flush=True)
            continue
        articles.append(article)
        labels.append('gwas')

    articles, labels = shuffle_articles_labels(articles, labels)

    # articles = articles[:100]
    # labels = labels[:100]

    print('Have ', len(articles), flush=True)
    counters = defaultdict(int)
    for label in labels:
        counters[label] += 1

    print('Counts per label:', counters, flush=True)
    print('Done converting', flush=True)
    classifier = text_classification.create_model(articles, labels, cross_val=cross_val, l1=l1)

    if save:
        print("SAVING TEXT FIELD RELEVANCE CLASSIFIER", flush=True)
        with open(out_dir + '/text_field_{}.pkl'.format(field_name), 'wb') as out_file:
            pickle.dump(classifier, out_file)
    else:
        print("NOT SAVING TEXT FIELD RELEVANCE CLASSIFIER", flush=True)


def train_topical_gene_classifier(out_dir, process_dir, allowed_pmids, forbidden_pmids,
                                  save, cross_val, l1,
                                  replace_phenos_with_nothing):
    with open(out_dir + '/dataset_meta.json') as file:
        dataset_info = json.load(file)

    with open(out_dir + '/omim.json') as file:
        omim = json.load(file)

    positives = set(dataset_info['positive_pmids']) - forbidden_pmids
    print('Total positives: ', len(positives), flush=True)

    articles = []
    labels = []

    for pmid, processed_article in convert_all_to_articles(process_dir, positives):
        if allowed_pmids is not None and pmid not in allowed_pmids:
            print("Skipping PMID %s because not in allowed PMIDs" % pmid, flush=True)
            continue
        omim_data = omim[str(pmid)]
        topical_genes = set(omim_data['topical_genes'])
        eids_in_document = set(processed_article['text']['gene_mentions'].keys())
        if len(topical_genes) > 3:
            continue

        for eid in eids_in_document:
            articles.append((processed_article, eid))
            labels.append(eid in topical_genes)

    articles, labels = shuffle_articles_labels(articles, labels)

    print('Have ', len(articles), flush=True)
    counters = defaultdict(int)

    for label in labels:
        counters[label] += 1

    print('Counts per label:', counters, flush=True)
    print('Done converting', flush=True)
    classifier = topical_gene.create_model(articles, labels, cross_val=cross_val, l1=l1)

    if save:
        print("SAVING TOPICAL GENE CLASSIFIER", flush=True)
        with open(out_dir + '/topical_gene.pkl', 'wb') as out_file:
            pickle.dump(classifier, out_file)
    else:
        print("NOT SAVING TOPICAL GENE CLASSIFIER", flush=True)


def train_variant_classifier(out_dir, process_dir, allowed_pmids, forbidden_pmids,
                             save, cross_val, replace_phenos_with_nothing):
    with open(out_dir + '/dataset_meta.json') as file:
        dataset_info = json.load(file)
    with open(out_dir + '/topical_gene.pkl', 'rb') as file:
        topical_gene_classifier = pickle.load(file)

    clinvar_db = clinvar.Clinvar()
    if allowed_pmids is not None:
        clinvar_pmids = list(set(dataset_info['clinvar_pmids']).intersection(allowed_pmids))
        print("Skipping %d PMIDs because not in allowed PMIDs" % (len(dataset_info['clinvar_pmids'])
                                                                  - len(clinvar_pmids)), flush=True)
    else:
        clinvar_pmids = dataset_info['clinvar_pmids']
    clinvar_pmids = set(clinvar_pmids) - forbidden_pmids

    print('Number of positive articles in the ClinVar PMIDs: ', len(clinvar_pmids), flush=True)

    featurized_variants = []
    labels = []

    for label, featurized_variant in \
            convert_all_to_variants(process_dir, clinvar_pmids,
                                    clinvar_db, topical_gene_classifier):
        labels.append(label)
        featurized_variants.append(featurized_variant)

    print('Have %d featurized variants' % len(featurized_variants), flush=True)

    featurized_variants, labels = shuffle_articles_labels(featurized_variants, labels)
    counters = defaultdict(int)
    for label in labels:
        counters[label] += 1

    print('Counts per label:', counters, flush=True)
    print('Done converting', flush=True)

    classifier = variant_classification.create_model(featurized_variants, labels, cross_val=cross_val)

    if save:
        print("SAVING VARIANT CLASSIFIER", flush=True)
        with open(out_dir + '/variant_classifier.pkl', 'wb') as out_file:
            pickle.dump(classifier, out_file)
    else:
        print("NOT SAVING VARIANT CLASSIFIER", flush=True)
