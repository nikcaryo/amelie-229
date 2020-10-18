import argparse

import multiprocessing

from contextlib import redirect_stdout

from . import parasol, pheno_extractor
from . import process_main

import json
import pickle
import os
import shutil
import sys
import random

import sqlite3

p_extractor = None


num_process_threads = multiprocessing.Value('i', 2)


def process_papers(args, pubmed_queue, processed_queue):
    extractors = {
        'p_extractor': p_extractor,
    }

    while True:
        pmid = pubmed_queue.get()
        if pmid is None:
            pubmed_queue.put(None)
            break
        else:
            process_status = process_main.process_paper(paper_dir=args.paper_dir,
                                                        process_dir=args.process_dir_xref_phenos,
                                                        extractors=extractors,
                                                        pubmed_info=None,
                                                        pmid=pmid,
                                                        fulltext_only=True)
            processed_queue.put((pmid, process_status))

    with num_process_threads.get_lock():
        num_process_threads.value -= 1

        if num_process_threads.value == 0:
            print('All processing work is done')
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
    return all(type(i)==str for i in fields)


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


def extract_xref_phenos_program():
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

    parser.add_argument('process_dir_xref_phenos',
        type=str,
        help='The directory where the processed papers are stored')

    parser.add_argument('pmids', type=str, nargs='+')

    args = parser.parse_args()
    with open(args.out_file, 'w') as out_file:
        print(args, file=out_file, flush=True)
        with redirect_stdout(out_file):
            if not os.path.isdir(args.out_dir):
                print("Out dir %s does not exist" % args.out_dir, file=out_file, flush=True)
                sys.exit(1)

            if not os.path.isdir(args.process_dir_xref_phenos):
                print("Paper dir %s does not exist" % args.process_dir_xref_phenos, file=out_file, flush=True)
                sys.exit(1)

            if not args.process_dir_xref_phenos.endswith('_xref_phenos'):
                print("As a cautionary measure, the process dir for "
                      "xref pheno extraction has to end with _xref_phenos")
                sys.exit(1)

            global p_extractor

            print("Loading pheno extractor", file=out_file, flush=True)
            p_extractor = pheno_extractor.PhenoExtractor(xref_synonyms=True)

            processes = []
            pubmed_queue = multiprocessing.Queue()
            processed_queue = multiprocessing.Queue()
            print('Num process:', num_process_threads.value, file=out_file, flush=True)

            for _ in range(num_process_threads.value):
                p = multiprocessing.Process(target=process_papers, args=(args, pubmed_queue, processed_queue))
                p.start()
                processes.append(p)

            print('Filling the queue', file=out_file, flush=True)
            for pmid in args.pmids:
                pubmed_queue.put(pmid)
            pubmed_queue.put(None)

            def extract(pmid, process_status):
                result = {}
                if not process_status['is_valid']:
                    return {'is_valid': False, 'reason': 'Process failure', 'process_status': process_status}
                with open(args.process_dir_xref_phenos + '/' + str(pmid) + '.pkl', 'rb') as file:
                    processed_article = pickle.load(file)
                result['phenotypes'] = list({x.hpo_id for processed_data in processed_article.values()
                                             for x in processed_data['phenotype_mentions']})
                result['is_valid'] = True
                return result

            while True:
                item = processed_queue.get()
                if item is None:
                    break
                else:
                    pmid, process_status = item
                    print('Processing ', pmid, file=out_file, flush=True)
                    result = extract(pmid, process_status)
                    print('Done processing ', pmid, file=out_file, flush=True)

                    extract_dir = args.out_dir + '/extracted_xref_phenos/'
                    os.makedirs(extract_dir, exist_ok=True)
                    with open(os.path.join(extract_dir, '%s.json' % pmid), 'w') as file:
                        json.dump(result, file)

            print('Joining processes', file=out_file, flush=True)
            for p in processes:
                p.join()
            print('Done processes', file=out_file, flush=True)


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def extract_all_xref_phenos_program():
    parser = argparse.ArgumentParser()

    parser.add_argument('amelie_knowledgebase',
                        type=str)

    parser.add_argument('out_dir',
        type=str)

    parser.add_argument('paper_dir',
        type=str)

    parser.add_argument('process_dir_xref_phenos',
        type=str)

    parser.add_argument('--conda_path', 
        type=str, 
        default='/cluster/u/ethanid/miniconda3/bin')

    parser.add_argument('--max_jobs',
                        type=int,
                        default=150)

    parser.add_argument('--delete_parasol_dir',
                        action="store_true",
                        default=False)

    parser.add_argument('--take_random',
                        type=int,
                        default=None)

    args = parser.parse_args()
    print(args)

    if not os.path.isdir(args.out_dir):
        print("Out dir %s does not exist" % args.out_dir)
        sys.exit(1)

    if not os.path.isdir(args.process_dir_xref_phenos):
        print("Process dir %s does not exist" % args.process_dir_xref_phenos)
        sys.exit(1)

    if not args.process_dir_xref_phenos.endswith('_xref_phenos'):
        print("As a cautionary measure, the process dir for "
              "xref pheno extraction has to end with _xref_phenos")
        sys.exit(1)

    if not os.path.isdir(args.paper_dir):
        print("Paper dir %s does not exist" % args.paper_dir)
        sys.exit(1)

    connection = sqlite3.connect(args.amelie_knowledgebase)
    cur = connection.cursor()
    qr = cur.execute("SELECT DISTINCT pmid FROM amelie_paper")
    all_pmids = [x[0] for x in qr.fetchall()]
    pmidss = chunks(all_pmids, 50)

    input_texts = [
        args.out_dir + ' ' +
        args.paper_dir + ' ' +
        args.process_dir_xref_phenos + ' ' +
        ' '.join(pmids)
        for pmids in pmidss
    ]

    print(input_texts[:3])

    random.shuffle(input_texts)
    if args.take_random is not None:
        input_texts = input_texts[:args.take_random]

    max_jobs = args.max_jobs
    num_cpu = 2
    num_mem = 4
    # intentionally not localized currently
    parasol_result = parasol.submit_jobs(args.conda_path, 'extract_xref_phenos', input_texts,
                                         max_jobs=max_jobs, num_cpu=num_cpu, num_mem=num_mem)

    if args.delete_parasol_dir:
        shutil.rmtree(parasol_result['run_directory'])
