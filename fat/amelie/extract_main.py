import argparse

import multiprocessing

from contextlib import redirect_stdout

from . import pubmed_batch, parasol, pheno_extractor, gene_extractor, variant_extractor, local_executor
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

num_find_relevant_papers = multiprocessing.Value('i', 8)
num_download_threads = multiprocessing.Value('i', 1)
num_process_threads = multiprocessing.Value('i', 8)


# WARNING !!! WARNING !!! WARNING:
#
# WHEN FIXING THIS FILE WITH IMPORTANT STUFF, MAKE SURE TO ALSO FIX
# extract_main_single_thread.py, WHICH IS USED ON THE AMELIE SERVER
# FOR CONTINUOUS UPDATES
###################################################################


pub = None
p_extractor = None
g_extractor = None


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
            pmid_queue.put(item)
            break
        else:
            for pmid, pubmed_info in item:
                if args.download_only:
                    if os.path.isfile(os.path.join(args.paper_dir, '%s.pdf' % pmid)):
                        print("Already downloaded %s, skipping" % pmid, flush=True)
                        continue
                text_results = {
                    'abstract': process_main.process_part_of_paper('abstract', pubmed_info['abstract'], extractors, pmid),
                    'title': process_main.process_part_of_paper('title', pubmed_info['title'], extractors, pmid),
                }
                
                features = text_classification.\
                    convert_to_text(text_results, use_main_text=False,
                                    replace_phenos_with_nothing=args.replace_phenos_with_nothing)
                probas = pubmed_only_relevance.predict_proba([features])
                assert probas.shape[0] == 1, probas.shape
                score = probas[0, 1]

                print("Tested title/abstract of PMID %s for relevance, result: %s" % (pmid, str(score > 0.5)), flush=True)
                if score > 0.5:
                    print("After testing for relevance, PMID\t%s\t%s\t%s" %
                          (pmid, pubmed_info['title'], pubmed_info['journal']), flush=True)
                    q.put((pmid, pubmed_info))
                
    with num_find_relevant_papers.get_lock():
        num_find_relevant_papers.value -= 1

        if num_find_relevant_papers.value == 0:
            print('All find relevant work is done', flush=True)
            q.put(None)


def download_papers(args, find_relevant_queue, downloaded_queue):
    config = {}

    if args.elsevier_key is not None:
        config['elsevierApiKey'] = args.elsevier_key

    if args.sfx_server is not None:
        config['crawlSfxServer'] = args.sfx_server

    while True:
        item = find_relevant_queue.get()
        if item is None:
            find_relevant_queue.put(None)
            break
        else:
            pmid, pubmed_info = item
            download_status = download_main.download_paper(args.paper_dir, config, pmid)
            downloaded_queue.put((pmid, pubmed_info, download_status))

    with num_download_threads.get_lock():
        num_download_threads.value -= 1

        if num_download_threads.value == 0:
            print('All downloading work is done', flush=True)
            downloaded_queue.put(None)


def dummy_download_papers_ta_only(args, find_relevant_queue, downloaded_queue):
    while True:
        item = find_relevant_queue.get()
        if item is None:
            find_relevant_queue.put(None)
            break
        else:
            pmid, pubmed_info = item
            download_status = download_main.dummy_download_paper_ta_only()
            downloaded_queue.put((pmid, pubmed_info, download_status))

    with num_download_threads.get_lock():
        num_download_threads.value -= 1

        if num_download_threads.value == 0:
            print('All downloading work is done', flush=True)
            downloaded_queue.put(None)


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
                downloaded_queue.put(None)
                break
            else:
                pmid, pubmed_info, download_status = item
                process_status = process_main.process_paper(args.paper_dir, args.process_dir, extractors,
                                                            pubmed_info, pmid)
                processed_queue.put((pmid, pubmed_info, download_status, process_status))

        with num_process_threads.get_lock():
            num_process_threads.value -= 1

            if num_process_threads.value == 0:
                print('All processing work is done', flush=True)
                processed_queue.put(None)


def process_papers_ta_only(args, downloaded_queue, processed_queue):
    extractors = {
        'p_extractor': p_extractor,
        'g_extractor': g_extractor,
    }

    while True:
        item = downloaded_queue.get()
        if item is None:
            downloaded_queue.put(None)
            break
        else:
            pmid, pubmed_info, download_status = item
            print("Got pmid %s from downloaded queue" % pmid, flush=True)
            process_status = process_main.process_title_abstract_only(args.paper_dir, args.process_dir, extractors,
                                                                      pubmed_info, pmid)
            processed_queue_item = (pmid, pubmed_info, download_status, process_status)
            print("Putting pmid %s into processed queue" % pmid, flush=True)
            processed_queue.put(processed_queue_item)

    with num_process_threads.get_lock():
        num_process_threads.value -= 1

        if num_process_threads.value == 0:
            print('All processing work is done', flush=True)
            processed_queue.put(None)


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
    parser.add_argument('--download_only',
                        action="store_true",
                        default=False)
    parser.add_argument('--title_abstract_relevance_classify_only',
                        action="store_true",
                        default=False)
    parser.add_argument('--elsevier_key', type=str, help='The pubmed id to download')
    parser.add_argument('--sfx_server',
                        type=str,
                        default='http://sul-sfx.stanford.edu/sfxlcl41')
    parser.add_argument('--replace_phenos_with_nothing', action='store_true', default=False)
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
                                                        replace_phenos_with_nothing=args.replace_phenos_with_nothing)
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
            pubmed_queue = multiprocessing.Queue()
            batch_items = list(batch.items())
            print('Filling the queue', file=out_file, flush=True)
            num_parts = 10
            items_per_part = (len(batch_items) + num_parts - 1) // num_parts
            batch_items = batch_items[args.part * items_per_part: (args.part + 1) * items_per_part]

            extract_dir = os.path.join(args.out_dir, args.extracted_subdir)

            for i in range(0, len(batch_items), 100):
                pubmed_queue.put(batch_items[i:i+100])

            pubmed_queue.put(None)

            processes = []

            global pub, p_extractor, g_extractor

            print("Loading pubtator", file=out_file, flush=True)
            pub = pubtator.Pubtator(args.out_dir + '/pubtator.db')
            print("Loading pheno extractor", file=out_file, flush=True)
            p_extractor = pheno_extractor.PhenoExtractor()
            print("Loading gene extractor", file=out_file, flush=True)
            g_extractor = gene_extractor.GeneExtractor(pub)

            find_relevant_queue = multiprocessing.Queue()
            downloaded_queue = multiprocessing.Queue()
            processed_queue = multiprocessing.Queue()

            print('Num find:', num_find_relevant_papers.value, file=out_file, flush=True)
            print('Num download:', num_download_threads.value, file=out_file, flush=True)
            print('Num process:', num_process_threads.value, file=out_file, flush=True)

            for _ in range(num_find_relevant_papers.value):
                p = multiprocessing.Process(target=find_relevant_papers, args=(args, pubmed_queue, find_relevant_queue))
                p.start()
                processes.append(p)
            create_download_processes(args, downloaded_queue, find_relevant_queue, processes)
            create_process_processes(args, downloaded_queue, processed_queue, processes)

            if args.title_abstract_relevance_classify_only:
                while True:
                    item = find_relevant_queue.get()
                    if item is None:
                        find_relevant_queue.put(None)
                        break
                return

            if args.download_only:
                while True:
                    item = downloaded_queue.get()
                    if item is None:
                        downloaded_queue.put(None)
                        break
                    else:
                        pmid, pubmed_info, downloaded_status = item
                        print("Relevant PMID %s, downloaded status: %s" % (pmid, downloaded_status),
                              file=out_file,
                              flush=True)
                return

            final_results = extract_from_processed_queue(args, batch_items, out_file, processed_queue, processes,
                                                         pubmed_queue)

            print('Joining processes', file=out_file, flush=True)
            for p in processes:
                p.join()
            print('Done processes', file=out_file, flush=True)

            os.makedirs(extract_dir, exist_ok=True)
            with open(os.path.join(extract_dir,
                                   str(args.update_year) + '-' +
                                   str(args.update_number) + '-' +
                                   str(args.part) + '.json'), 'w') as file:
                json.dump(final_results, file)


def extract_from_processed_queue(args, batch_items, out_file, processed_queue, processes, pubmed_queue):
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


def create_process_processes(args, downloaded_queue, processed_queue, processes):
    if not args.download_only and not args.title_abstract_relevance_classify_only:
        if args.title_abstract_only:
            for _ in range(num_process_threads.value):
                p = multiprocessing.Process(target=process_papers_ta_only,
                                            args=(args, downloaded_queue, processed_queue))
                p.start()
                processes.append(p)
        else:
            for _ in range(num_process_threads.value):
                p = multiprocessing.Process(target=process_papers, args=(args, downloaded_queue, processed_queue))
                p.start()
                processes.append(p)


def create_download_processes(args, downloaded_queue, find_relevant_queue, processes):
    if not args.title_abstract_relevance_classify_only:
        if args.title_abstract_only:
            for _ in range(num_download_threads.value):
                p = multiprocessing.Process(target=dummy_download_papers_ta_only, args=(args,
                                                                                        find_relevant_queue,
                                                                                        downloaded_queue))
                p.start()
                processes.append(p)
        else:
            for _ in range(num_download_threads.value):
                p = multiprocessing.Process(target=download_papers, args=(args, find_relevant_queue, downloaded_queue))
                p.start()
                processes.append(p)


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
        ("--download_only " if args.download_only else " ") +
        ("--title_abstract_relevance_classify_only " if args.title_abstract_relevance_classify_only else " ") +
        ("--replace_phenos_with_nothing " if args.replace_phenos_with_nothing else " ") +
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

    if args.download_only:
        max_jobs = args.max_jobs
        num_cpu = 8
        num_mem = 20
    else:
        max_jobs = args.max_jobs
        num_cpu = 8
        num_mem = 20

    print("Have %d updates*parts to process" % len(input_texts), flush=True)
    if args.local_processing:
        if len(input_texts) > 1000:
            print("Over 1000 updates (%d) to process, do not run this locally" %  len(input_texts), flush=True)
            sys.exit(1)
        parasol_result = local_executor.submit_jobs(args.conda_path, 'extract_papers', input_texts)
    else:
        parasol_result = parasol.submit_jobs(args.conda_path, 'extract_papers', input_texts,
                                             max_jobs=max_jobs, num_cpu=num_cpu, num_mem=num_mem)

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
    parser.add_argument('--download_only',
                        action="store_true",
                        default=False)
    parser.add_argument('--title_abstract_relevance_classify_only',
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
    parser.add_argument('--replace_phenos_with_nothing', action='store_true', default=False)
    parser.add_argument('--local_processing', action='store_true', default=False)
    parser.add_argument('--last_output_tmpdir', action='store_true', default=False)
    args = parser.parse_args()
    return args
