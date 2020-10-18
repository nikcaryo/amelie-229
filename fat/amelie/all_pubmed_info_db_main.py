import argparse

from contextlib import redirect_stdout

import os
import shutil
import sys
import random
import pickle

from . import pubmed_batch, parasol, pheno_extractor, gene_extractor, local_executor, process_main, text_classification
from .dataloaders import pubtator


pub = None
p_extractor = None
g_extractor = None


def is_paper_relevant(args, pmid, paper_info, pubmed_only_relevance_clf):
    extractors = {
        'p_extractor': p_extractor,
        'g_extractor': g_extractor
    }

    if 'abstract' not in paper_info:
        print(paper_info)
        print(paper_info, file=sys.stderr)

    text_results = {
        'abstract': process_main.process_part_of_paper('abstract', paper_info['abstract'], extractors, pmid),
        'title': process_main.process_part_of_paper('title', paper_info['title'], extractors, pmid),
    }

    features = text_classification. \
        convert_to_text(text_results, use_main_text=False,
                        replace_phenos_with_nothing=False)
    probas = pubmed_only_relevance_clf.predict_proba([features])
    assert probas.shape[0] == 1, probas.shape
    score = probas[0, 1]
    return score > 0.5


def paper_downloaded(args, pmid):
    input_location = os.path.join(args.paper_dir, '%s.pdf' % pmid)
    return os.path.isfile(input_location) and os.stat(input_location).st_size > 0


def pubmed_info_db_program():
    parser = argparse.ArgumentParser()
    parser.add_argument('out_file', type=str, help='The file to output debug information to')
    parser.add_argument('pubmed_info_dir', type=str)
    parser.add_argument('update_number', type=int, help='The pubmed update id')
    parser.add_argument('--classify_ta_relevance', action="store_true", default=False)
    parser.add_argument('--lookup_pdf_exists', action="store_true", default=False)
    parser.add_argument('--paper_dir', type=str, default=None)
    parser.add_argument('--out_dir', type=str, default=None)
    args = parser.parse_args()

    with open(args.out_file, 'w') as out_file:
        print(args, file=out_file, flush=True)

        with redirect_stdout(out_file):

            if args.classify_ta_relevance:
                global pub, p_extractor, g_extractor

                print("Loading pubtator", file=out_file, flush=True)
                pub = pubtator.Pubtator(args.out_dir + '/pubtator.db')
                print("Loading pheno extractor", file=out_file, flush=True)
                p_extractor = pheno_extractor.PhenoExtractor()
                print("Loading gene extractor", file=out_file, flush=True)
                g_extractor = gene_extractor.GeneExtractor(pub)
                with open(args.out_dir + '/pubmed_only_relevance.pkl', 'rb') as file:
                    pubmed_only_relevance_clf = pickle.load(file)
            else:
                pubmed_only_relevance_clf = None

            if not os.path.isdir(args.pubmed_info_dir):
                print("Pubmed info dir %s does not exist" % args.pubmed_info_dir, file=out_file, flush=True)
                sys.exit(1)

            print("Connecting to pubmed to find out latest numbers", file=out_file, flush=True)
            baseline_num_pieces, year = pubmed_batch.get_num_pieces_and_year(pubmed_type='baseline')
            updates_num_pieces, year = pubmed_batch.get_num_pieces_and_year(pubmed_type='updatefiles')
            print("Got latest numbers from PubMed")

            if 1 <= args.update_number <= baseline_num_pieces:
                batch = pubmed_batch.process_pubmed_batch_file(args.update_number, year, pubmed_type='baseline')
            elif (baseline_num_pieces + 1) <= args.update_number <= updates_num_pieces:
                batch = pubmed_batch.process_pubmed_batch_file(args.update_number, year, pubmed_type='updatefiles')
            else:
                print('Provided batch number,', args.update_number, ' is not a valid batch number', file=out_file, flush=True)
                sys.exit(1)
                return  # to stop pycharm from complaining about not initialized batch variable

            def extract_some_info(pubmed_info):
                rv = []
                for pmid, info in pubmed_info.items():
                    if args.classify_ta_relevance:
                        if is_paper_relevant(args=args, pmid=pmid,
                                             paper_info=info,
                                             pubmed_only_relevance_clf=pubmed_only_relevance_clf):
                            ta_relevance = "Relevant"
                        else:
                            ta_relevance = "Irrelevant"
                    else:
                        ta_relevance = "Not classified"
                    if args.lookup_pdf_exists:
                        if paper_downloaded(args=args, pmid=pmid):
                            pdf_exists = "Downloaded"
                        else:
                            pdf_exists = "Not downloaded"
                    else:
                        pdf_exists = "Not checked"
                    rv.append([pmid,
                               info['title'].replace('\t', ' ').replace('\n', ' '),
                               info['publish_date'].replace('\t', ' ').replace('\n', ' ') if info['publish_date'] is not None else "None",
                               info['journal'].replace('\t', ' ').replace('\n', ' '),
                               '|^|'.join(info['authors']).replace('\t', ' ').replace('\n', ' '),
                               ta_relevance, pdf_exists])
                return rv

            os.makedirs(args.pubmed_info_dir, exist_ok=True)
            elems = extract_some_info(batch)
            with open(os.path.join(args.pubmed_info_dir, str(args.update_number) + '.tsv'), 'w') as file:
                for elem in elems:
                    print('\t'.join(elem), file=file)


def all_pubmed_info_db_program():
    parser = argparse.ArgumentParser()
    parser.add_argument('pubmed_info_dir', type=str)
    parser.add_argument('--conda_path', type=str, default='/cluster/u/jbirgmei/miniconda3/bin')
    parser.add_argument('--max_jobs', type=int, default=200)
    parser.add_argument('--local_processing', action='store_true', default=False)
    parser.add_argument('--last_output_tmpdir', action='store_true', default=False)
    parser.add_argument('--delete_parasol_dir', action="store_true", default=False)
    # these options make the thing insanely slow and memory heavy
    parser.add_argument('--classify_ta_relevance', action="store_true", default=False)
    parser.add_argument('--lookup_pdf_exists', action="store_true", default=False)
    parser.add_argument('--paper_dir', type=str, default=None)
    parser.add_argument('--out_dir', type=str, default=None)

    args = parser.parse_args()

    assert not args.classify_ta_relevance or (args.paper_dir and args.out_dir), args

    print(args, flush=True)

    updates_num_pieces, year = pubmed_batch.get_num_pieces_and_year(pubmed_type='baseline')
    update_numbers = []

    for update_number in range(1, updates_num_pieces + 1):
        extract_dir = args.pubmed_info_dir
        extract_file = os.path.join(extract_dir, str(update_number) + '.tsv')
        if not os.path.isfile(extract_file):
            update_numbers.append(update_number)

    input_texts = [
        args.pubmed_info_dir + ' ' +
        str(update_number) + ' ' +
        ('--classify_ta_relevance ' if args.classify_ta_relevance else "") +
        ('--lookup_pdf_exists ' if args.lookup_pdf_exists else "") +
        ('--paper_dir ' + args.paper_dir + ' ' if args.paper_dir else "") +
        ('--out_dir ' + args.out_dir + ' ' if args.out_dir else "")
        for update_number in update_numbers
    ]

    print(input_texts[:3], flush=True)

    random.shuffle(input_texts)

    max_jobs = args.max_jobs
    num_cpu = 1
    if args.classify_ta_relevance:
        num_mem = 8
    else:
        num_mem = 4

    if args.local_processing:
        parasol_result = local_executor.submit_jobs(args.conda_path, 'pubmed_info_db', input_texts)
    else:
        parasol_result = parasol.submit_jobs(args.conda_path, 'pubmed_info_db', input_texts,
                                             max_jobs=max_jobs, num_cpu=num_cpu, num_mem=num_mem)

    if args.delete_parasol_dir:
        shutil.rmtree(parasol_result['run_directory'])

    if args.last_output_tmpdir:
        print(flush=True)
        print(parasol_result['run_directory'], flush=True)
