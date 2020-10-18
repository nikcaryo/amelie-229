import argparse

from contextlib import redirect_stdout

from . import pubmed_batch, pheno_extractor, gene_extractor, variant_extractor, local_executor
from .refseq_extractor import RefseqExtractor
from . import position_extractor, process_main, text_classification, download_main, variant_classification
from .dataloaders import pubtator

import json
import pickle
import os
import shutil
import sys
import random

import urllib.error

pub = None
p_extractor = None
g_extractor = None


class Queue:

    def __init__(self):
        self.queue = []

    def put(self, elem):
        self.queue.append(elem)

    def get(self):
        if len(self.queue) > 0:
            return self.queue.pop(0)
        return None

    def qsize(self):
        return len(self.queue)


def find_relevant_papers(args, pmid_queue, q):
    with open(args.out_dir + '/pubmed_only_relevance.pkl', 'rb') as file:
        pubmed_only_relevance = pickle.load(file)

    extractors = {
        'p_extractor': p_extractor,
        'g_extractor': g_extractor
    }

    while True:
        item = pmid_queue.get()
        if item is None:
            break
        else:
            for pmid, pubmed_info in item:
                text_results = {
                    'abstract': process_main.process_part_of_paper('abstract', pubmed_info['abstract'], extractors, pmid),
                    'title': process_main.process_part_of_paper('title', pubmed_info['title'], extractors, pmid),
                }

                features = text_classification.\
                    convert_to_text(text_results, use_main_text=False,
                                    replace_phenos_with_nothing=False)
                probas = pubmed_only_relevance.predict_proba([features])
                assert probas.shape[0] == 1, probas.shape
                score = probas[0, 1]

                print("Tested title/abstract of PMID %s for relevance, result: %s" % (pmid, str(score > 0.5)), flush=True)
                if score > 0.5:
                    print("After testing for relevance, PMID\t%s\t%s\t%s" %
                          (pmid, pubmed_info['title'], pubmed_info['journal']), flush=True)
                    q.put((pmid, pubmed_info))


def download_papers(args, find_relevant_queue, downloaded_queue):
    config = {}

    if args.elsevier_key is not None:
        config['elsevierApiKey'] = args.elsevier_key

    if args.sfx_server is not None:
        config['crawlSfxServer'] = args.sfx_server

    while True:
        item = find_relevant_queue.get()
        if item is None:
            break
        else:
            pmid, pubmed_info = item
            download_status = download_main.download_paper(args.paper_dir, config, pmid)
            downloaded_queue.put((pmid, pubmed_info, download_status))


def process_papers(args, downloaded_queue, processed_queue):
    v_extractor = variant_extractor.VariantExtractor(args.pubmunch_data_dir)
    refseq_extractor = RefseqExtractor()
    with position_extractor.PositionExtractor() as pos_extractor:
        extractors = {
            'p_extractor': p_extractor,
            'g_extractor': g_extractor,
            'v_extractor': v_extractor,
            'pos_extractor': pos_extractor,
            'refseq_extractor': refseq_extractor,
        }

        while True:
            item = downloaded_queue.get()
            if item is None:
                break
            else:
                pmid, pubmed_info, download_status = item
                process_status = process_main.process_paper(args.paper_dir, args.process_dir, extractors,
                                                            pubmed_info, pmid)
                processed_queue.put((pmid, pubmed_info, download_status, process_status))


# From https://stackoverflow.com/questions/33181170/how-to-convert-a-nested-namedtuple-to-a-dict
def isnamedtupleinstance(x):
    _type = type(x)
    bases = _type.__bases__
    if len(bases) != 1 or bases[0] != tuple:
        return False
    fields = getattr(_type, '_fields', None)
    if not isinstance(fields, tuple):
        return False
    return all(type(i) == str for i in fields)


def unpack(obj):
    if isinstance(obj, dict):
        return {key: unpack(value) for key, value in obj.items()}
    elif obj is None:
        return None
    elif isinstance(obj, list):
        return [unpack(value) for value in obj]
    elif isnamedtupleinstance(obj):
        return {key: unpack(value) for key, value in obj._asdict().items()}
    elif isinstance(obj, tuple):
        return tuple(unpack(value) for value in obj)
    elif isinstance(obj, str):
        return obj
    elif isinstance(obj, float):
        return obj
    elif isinstance(obj, int):
        return obj
    else:
        return {s: getattr(obj, s, None) for s in obj.__slots__}


def parse_extract_papers_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('out_file',
                        type=str,
                        help='The file to output debug information to')
    parser.add_argument('out_dir',
                        type=str,
                        help='The outdirectory where models and extracts are stored')
    parser.add_argument('paper_dir',
                        type=str,
                        help='The directory where the downloaded papers are stored')
    parser.add_argument('process_dir',
                        type=str,
                        help='The directory where the processed papers are stored')
    parser.add_argument('update_year', type=str)
    parser.add_argument('update_number', type=int, help='The pubmed update id')
    parser.add_argument('part', type=int,
                        help='The pubmed updates are split into 10 parts. This is the part to process.'
                        )
    parser.add_argument('pubmed_type', type=str, help='"baseline" or "updatefiles"')
    parser.add_argument('--title_abstract_only',
                        action="store_true",
                        default=False)
    parser.add_argument('extracted_subdir', type=str)
    parser.add_argument('--pubmunch_data_dir',
                        type=str,
                        default='/cluster/u/jbirgmei/shallowfloat/amelie_ethan/pubmunchData',
                        help='The location of your pubmunch data directory'
                        )
    parser.add_argument('--elsevier_key', type=str, help='The pubmed id to download')
    parser.add_argument('--sfx_server',
                        type=str,
                        default='http://sul-sfx.stanford.edu/sfxlcl41')
    args = parser.parse_args()
    return args


def extract_ml_scores(pmid, pubmed_info, download_status, process_status,
                      args,
                      pdf_relevance_classifier, text_field_inheritance_modes_classifier,
                      text_field_variant_types_classifier, topical_gene_classifier,
                      variant_classifier):
    result = dict()
    result['publish_date'] = pubmed_info['publish_date']
    result['journal'] = pubmed_info['journal']
    result['title'] = pubmed_info['title']
    result['authors'] = ', '.join(pubmed_info['authors'])

    if not download_status['is_valid']:
        return {'is_valid': False, 'reason': 'Download failure', 'download_status': download_status}

    if not process_status['is_valid']:
        return {'is_valid': False, 'reason': 'Process failure', 'process_status': process_status}

    with open(args.process_dir + '/' + str(pmid) + '.pkl', 'rb') as file:
        processed_article = pickle.load(file)
    text_features = text_classification.convert_to_text(processed_article, use_main_text=True,
                                                        replace_phenos_with_nothing=False)
    probas = pdf_relevance_classifier.predict_proba([text_features])
    assert probas.shape[0] == 1, probas.shape
    result['relevance'] = probas[0, 1]

    def process_text_field(name, classifier):
        probs = classifier.predict_proba([text_features])
        assert probs.shape[0] == 1, probs.shape
        probs = probs[0, :]
        result[name] = {
            class_name: probs[i] for i, class_name in enumerate(classifier.classes_)
        }

    process_text_field('inheritance_modes', text_field_inheritance_modes_classifier)
    process_text_field('variant_types', text_field_variant_types_classifier)

    result['topical_genes'] = {}
    for eid in processed_article['text']['gene_mentions']:
        probas = topical_gene_classifier.predict_proba([(processed_article, eid)])
        assert probas.shape[0] == 1, probas.shape
        result['topical_genes'][eid] = probas[0, 1]

    result['variants'] = []
    if len(processed_article['text']['variant_mentions']) > 0 and variant_classifier is not None:
        print("Number of variant "
              "mentions in PMID %s: %d" % (pmid,
                                           len(processed_article['text']['variant_mentions'])), flush=True)
        if len(processed_article['text']['variant_mentions']) > 1000:
            print("Have too many variant mentions in PMID %s, skipping" % pmid, flush=True)
        else:
            v_feature_extractor = variant_classification.DocumentVariantFeatureExtractor(processed_article,
                                                                                         topical_gene_classifier)
            for variant in processed_article['text']['variant_mentions']:
                variant_rows = v_feature_extractor.featurize(variant)
                if len(variant_rows) > 0:
                    probas = variant_classifier.predict_proba(variant_rows)
                    variant_val = max(probas[:, 1])
                    result['variants'].append({'variant': unpack(variant), 'variant_score': variant_val})
    result['phenotypes'] = list({x.hpo_id for processed_data in processed_article.values()
                                 for x in processed_data['phenotype_mentions']})
    result['is_valid'] = True
    return result


def extract_papers_program():
    args = parse_extract_papers_args()
    with open(args.out_file, 'w') as out_file:
        print(args, file=out_file, flush=True)
        with redirect_stdout(out_file):
            check_dirs_exist(args, out_file)

            batch = pubmed_batch.process_pubmed_batch_file(args.update_number,
                                                           year=args.update_year,
                                                           pubmed_type=args.pubmed_type)
            pubmed_queue = Queue()
            batch_items = list(batch.items())
            print('Filling the queue', file=out_file, flush=True)
            num_parts = 10
            items_per_part = (len(batch_items) + num_parts - 1) // num_parts
            batch_items = batch_items[args.part * items_per_part: (args.part + 1) * items_per_part]

            extract_dir = os.path.join(args.out_dir, args.extracted_subdir)

            for i in range(0, len(batch_items), 100):
                pubmed_queue.put(batch_items[i:i+100])

            pubmed_queue.put(None)

            global pub, p_extractor, g_extractor

            print("Loading pubtator", file=out_file, flush=True)
            pub = pubtator.Pubtator(args.out_dir + '/pubtator.db')
            print("Loading pheno extractor", file=out_file, flush=True)
            p_extractor = pheno_extractor.PhenoExtractor()
            print("Loading gene extractor", file=out_file, flush=True)
            g_extractor = gene_extractor.GeneExtractor(pub)

            find_relevant_queue = Queue()
            downloaded_queue = Queue()
            processed_queue = Queue()

            find_relevant_papers(args, pubmed_queue, find_relevant_queue)
            download_papers(args, find_relevant_queue, downloaded_queue)
            process_papers(args, downloaded_queue, processed_queue)
            final_results = extract_from_processed_queue(args, batch_items, out_file, processed_queue, pubmed_queue)

            print('Done with final results', file=out_file, flush=True)

            os.makedirs(extract_dir, exist_ok=True)
            with open(os.path.join(extract_dir,
                                   str(args.update_year) + '-' +
                                   str(args.update_number) + '-' +
                                   str(args.part) + '.json'), 'w') as file:
                json.dump(final_results, file)


def extract_from_processed_queue(args, batch_items, out_file, processed_queue, pubmed_queue):
    print('Setting up extractors', flush=True)
    with open(args.out_dir + '/pdf_relevance.pkl', 'rb') as file:
        pdf_relevance_classifier = pickle.load(file)
    with open(args.out_dir + '/text_field_inheritance_modes.pkl', 'rb') as file:
        text_field_inheritance_modes_classifier = pickle.load(file)
    with open(args.out_dir + '/text_field_variant_types.pkl', 'rb') as file:
        text_field_variant_types_classifier = pickle.load(file)
    with open(args.out_dir + '/topical_gene.pkl', 'rb') as file:
        topical_gene_classifier = pickle.load(file)
    if args.title_abstract_only:
        variant_classifier = None
    else:
        with open(args.out_dir + '/variant_classifier.pkl', 'rb') as file:
            variant_classifier = pickle.load(file)
    final_results = {}
    while True:
        item = processed_queue.get()
        if item is None:
            break
        else:
            pmid, pubmed_info, download_status, process_status = item
            print('Processing ', pmid, pubmed_queue.qsize(), len(batch_items) // 100, file=out_file, flush=True)
            result = extract_ml_scores(pmid, pubmed_info, download_status, process_status,
                                       args=args,
                                       pdf_relevance_classifier=pdf_relevance_classifier,
                                       text_field_inheritance_modes_classifier=text_field_inheritance_modes_classifier,
                                       text_field_variant_types_classifier=text_field_variant_types_classifier,
                                       topical_gene_classifier=topical_gene_classifier,
                                       variant_classifier=variant_classifier)
            print('Done processing ', pmid, file=out_file, flush=True)
            final_results[pmid] = result
    return final_results


def check_dirs_exist(args, out_file):
    if not os.path.isdir(args.out_dir):
        print("Out dir %s does not exist" % args.out_dir, file=out_file, flush=True)
        sys.exit(1)
    if not args.title_abstract_only:
        if not os.path.isdir(args.paper_dir):
            print("Paper dir %s does not exist" % args.paper_dir, file=out_file, flush=True)
            sys.exit(1)
    if not os.path.isdir(args.process_dir):
        print("Paper dir %s does not exist" % args.process_dir, file=out_file, flush=True)
        sys.exit(1)
    if not os.path.isdir(args.pubmunch_data_dir):
        print("Paper dir %s does not exist" % args.pubmunch_data_dir, file=out_file, flush=True)
        sys.exit(1)


def extract_all_papers_program():
    args = create_extract_all_papers_args()
    print(args, flush=True)

    check_extract_all_papers_dirs(args)

    if not args.no_pubtator_download:
        print("Downloading pubtator", flush=True)
        try:
            pubtator.download_pubtator(args.out_dir + '/pubtator.db')
            print("Download completed", flush=True)
        except urllib.error.URLError:
            print("PubTator download error", flush=True)


    baseline_num_pieces, baseline_year = pubmed_batch.get_num_pieces_and_year(pubmed_type='baseline')
    updates_num_pieces, update_year = pubmed_batch.get_num_pieces_and_year(pubmed_type='updatefiles')
    assert baseline_year == update_year, (baseline_year, update_year)

    parts = []

    if args.max_update_number is not None:
        updates_num_pieces = min(args.max_update_number, updates_num_pieces)

    for update_number in range(1, updates_num_pieces + 1):
        if 1 <= update_number <= baseline_num_pieces:
            pubmed_type = 'baseline'
        elif (baseline_num_pieces + 1) <= update_number <= updates_num_pieces:
            pubmed_type = 'updatefiles'
        else:
            pubmed_type = None
            print('Provided batch number ', args.update_number, ' is not a valid batch number', flush=True)
            sys.exit(1)
        for part in range(10):
            extract_dir = os.path.join(args.out_dir, args.extracted_subdir)
            extract_file = os.path.join(extract_dir,
                                        str(update_year) + '-' + str(update_number) + '-' + str(part) + '.json')
            if not os.path.isfile(extract_file):
                parts.append((update_number, part, pubmed_type))

    input_texts = [
        '--pubmunch_data_dir ' + args.pubmunch_data_dir + ' ' +
        '--elsevier_key ' + args.elsevier_key + ' ' +
        '--sfx_server ' + args.sfx_server + ' ' +
        args.out_dir + ' ' +
        args.paper_dir + ' ' +
        args.process_dir + ' ' +
        str(update_year) + ' ' +
        str(update_number) + ' ' +
        str(part) + ' ' +
        str(pubmed_type) + ' ' +
        args.extracted_subdir + ' '
        for (update_number, part, pubmed_type) in parts]

    if args.title_abstract_only:
        input_texts = [x + ' --title_abstract_only' for x in input_texts]

    print(input_texts[:3], flush=True)

    random.shuffle(input_texts)
    if args.take_random is not None:
        input_texts = input_texts[:args.take_random]

    print("Have %d updates*parts to process" % len(input_texts), flush=True)
    if len(input_texts) > 1000:
        print("Over 1000 updates (%d) to process, do not run this locally" % len(input_texts), flush=True)
        sys.exit(1)
    parasol_result = local_executor.submit_jobs(args.conda_path, 'extract_papers_single_thread', input_texts)

    if args.delete_parasol_dir:
        shutil.rmtree(parasol_result['run_directory'])

    if args.last_output_tmpdir:
        print(flush=True)
        print(parasol_result['run_directory'], flush=True)


def check_extract_all_papers_dirs(args):
    if not os.path.isdir(args.out_dir):
        print("Out dir %s does not exist" % args.out_dir, flush=True)
        sys.exit(1)
    if not os.path.isdir(args.process_dir):
        print("Process dir %s does not exist" % args.process_dir, flush=True)
        sys.exit(1)
    if not args.title_abstract_only:
        if not os.path.isdir(args.paper_dir):
            print("Paper dir %s does not exist" % args.paper_dir, flush=True)
            sys.exit(1)


def create_extract_all_papers_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('out_dir',
                        type=str,
                        help='The name of the database to output')
    parser.add_argument('paper_dir',
                        type=str)
    parser.add_argument('process_dir',
                        type=str,
                        default='processed_papers',
                        help='The name of the database to output')
    parser.add_argument('--pubmunch_data_dir',
                        type=str,
                        default='/cluster/u/ethanid/pubmunchData',
                        help='The location of your pubmed baseline dump'
                        )
    parser.add_argument('--conda_path',
                        type=str,
                        default='/cluster/u/ethanid/miniconda3/bin',
                        help='The location of your pubmed baseline dump')
    parser.add_argument('--elsevier_key',
                        type=str,
                        default='dummy_key',
                        help='The pubmed id to download')
    parser.add_argument('--max_jobs',
                        type=int,
                        default=50)
    parser.add_argument('--sfx_server',
                        type=str,
                        default='http://sul-sfx.stanford.edu/sfxlcl41')
    parser.add_argument('--no_pubtator_download',
                        action='store_true',
                        default=False)
    parser.add_argument('--delete_parasol_dir',
                        action="store_true",
                        default=False)
    parser.add_argument('--take_random',
                        type=int,
                        default=None)
    parser.add_argument('--title_abstract_only',
                        action='store_true',
                        default=False)
    parser.add_argument('--max_update_number',
                        type=int,
                        default=None)
    parser.add_argument('--extracted_subdir', type=str, default='extracted/')
    parser.add_argument('--last_output_tmpdir', action='store_true', default=False)
    args = parser.parse_args()
    return args
