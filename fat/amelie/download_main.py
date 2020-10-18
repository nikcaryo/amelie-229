import argparse
import json
import multiprocessing
import os
import random
import sys
import traceback

import constdb
import pubMunch

from amelie.dataloaders import clinvar
from .dataloaders import omim, pubtator, gwas
from . import parasol
from . import pubmed_batch


def exception_to_str(ex):
    return ''.join(traceback.format_exception(etype=type(ex), value=ex, tb=ex.__traceback__))


class DownloadPubmedFunc:
    def __init__(self, year):
        self.year = year

    def __call__(self, current_index):
        return current_index, pubmed_batch.process_pubmed_batch_file(current_index, self.year)


def download_pubmed_baseline(out_file, num_threads):
    num_pieces, year = pubmed_batch.get_num_pieces_and_year()

    print(year, num_pieces)

    pool = multiprocessing.Pool(num_threads)

    pubmed_files_to_pull = range(1, num_pieces + 1)

    with constdb.create(out_file) as output:

        processed_info = {}

        for i, (index, data) in enumerate(pool.imap_unordered(DownloadPubmedFunc(year), pubmed_files_to_pull)):
            print('Processed ', i, index)
            for pmid, article_data in data.items():
                if int(pmid) in processed_info:
                    print('Duplicate?', pmid, i, index, processed_info[int(pmid)])
                processed_info[int(pmid)] = index
                output.add(int(pmid), json.dumps(article_data).encode('utf8'))


class DownloadMimFunc:
    def __init__(self, api_key):
        self.api_key = api_key

    def __call__(self, name):
        return name, omim.get_omin_file(name, self.api_key)


def download_omim(out_file, omim_key, hgmd_location, num_threads):
    mim_map = omim.get_mim_names_to_ensemble()

    mim_names = list(mim_map.keys())

    pool = multiprocessing.Pool(num_threads)

    omim_map = {}

    for i, (name, data) in enumerate(pool.imap_unordered(DownloadMimFunc(omim_key), mim_names)):
        if i % 100 == 0:
            print('Processed {} {}/{}'.format(name, i, len(mim_names)))
        omim_map[name] = data

    with open(out_file, 'w') as outfile:
        json.dump(omim.get_disease_database(omim_map, hgmd_location, mim_map), outfile)


def download_data_program():
    parser = argparse.ArgumentParser(description='The amelie tool downloads and processes ')

    parser.add_argument('out_dir', type=str, 
        help='The name of the database to output')

    parser.add_argument('--omim_key', type=str, default='jDJzYZ0jQAuhiFF-c0s32g', 
        help='Your API key to access the OMIM database')

    parser.add_argument('--hgmd_location', 
        type=str, 
        default='/cluster/data/labResources/medgen/hgmd/HGMD_PRO_2018.1/allmut-table-2018.1.tsv'
    )

    parser.add_argument('--num_threads', type=int, default=10, 
        help='Your API key to access the OMIM database')

    args = parser.parse_args()
    print(args)

    os.mkdir(args.out_dir)

    download_omim(args.out_dir + '/omim.json', args.omim_key, args.hgmd_location, args.num_threads)
    pubtator.download_pubtator(args.out_dir + '/pubtator.db')
    download_pubmed_baseline(args.out_dir + '/pubmed.db', args.num_threads)


def download_paper(out_directory, config, pmid):
    pdf_output_location = out_directory + '/' + str(pmid) + '.pdf'

    if os.path.isfile(pdf_output_location):
        return {
            'is_valid': True,
            'status': 'Already downloaded'
        }

    with open(pdf_output_location, 'wb') as outfile:
        try:
            # return {
            #     'is_valid': False,
            #     'message': 'Downloads are disabled'
            # }

            pdf = pubMunch.download_pmid(pmid, config=config)

            if pdf is not None:
                    outfile.write(pdf)
                    return {
                        'is_valid': True,
                        'status': 'Downloaded'
                    }
            else:
                return {
                    'is_valid': False,
                    'message': 'Could not download'
                }
        except Exception as e:
            print(e)
            return {
                'is_valid': False,
                'status': 'Error',
                'error_message': exception_to_str(e)
            }


def dummy_download_paper_ta_only():
    return {
        'is_valid': True,
        'status': 'Title/abstract only classification'
    }


def download_papers_program():
    parser = argparse.ArgumentParser()

    parser.add_argument('status_file', 
        type=str, 
        help='')

    parser.add_argument('out_dir',
        type=str, 
        help='The name of the database to output')

    parser.add_argument('pubmed_ids', 
        type=int,
        nargs='*',
        help='The pubmed ids to download')

    parser.add_argument('--elsevier_key', type=str, help='The pubmed id to download')

    parser.add_argument('--sfx_server',
        type=str,
        default='http://sul-sfx.stanford.edu/sfxlcl41')

    args = parser.parse_args()

    if not os.path.isdir(args.out_dir):
        print("Out dir %s does not exist" % args.out_dir)
        sys.exit(1)

    print(args)

    config = {}

    if args.elsevier_key is not None:
        config['elsevierApiKey'] = args.elsevier_key

    if args.sfx_server is not None:
        config['crawlSfxServer'] = args.sfx_server

    results = {}

    for pmid in args.pubmed_ids:
        results[pmid] = download_paper(args.out_dir, config, pmid)

    with open(args.status_file, 'w') as outfile:
        json.dump(results, outfile)


def download_all_papers_program():
    parser = argparse.ArgumentParser()

    parser.add_argument('out_dir',
        type=str,
        help='The location where to store the running data'
    )

    parser.add_argument('paper_dir',
        type=str, 
        help='The name of the database to output')

    parser.add_argument('--elsevier_key', 
        type=str, 
        default='2aaccd20dfd0e61f862238096bb05cf6',
        help='The pubmed id to download')

    parser.add_argument('--sfx_server',
        type=str,
        default='http://sul-sfx.stanford.edu/sfxlcl41')

    parser.add_argument('--conda_path', 
        type=str, 
        default='/cluster/u/ethanid/miniconda3/bin', 
        help='The location of your pubmed baseline dump')

    parser.add_argument('--negative_pmids_list',
                        type=str,
                        default=None)

    args = parser.parse_args()
    print(args)

    if not os.path.isdir(args.out_dir):
        print("Out dir %s does not exist" % args.out_dir)
        sys.exit(1)

    if not os.path.isdir(args.paper_dir):
        print("Paper dir %s does not exist" % args.paper_dir)
        sys.exit(1)

    random.seed(4537756)
    
    with open(args.out_dir + '/omim.json') as omim_dump:
        omim = json.load(omim_dump)
        all_positive_pmids = {int(a) for a in omim.keys()}

    print('Num total positives ', len(all_positive_pmids))

    with constdb.read(args.out_dir + '/pubmed.db') as pubmed_dump:
        pubmed_pmids = set(pubmed_dump.keys())

    positive_pmids = pubmed_pmids & all_positive_pmids

    print('Total positives: ', len(positive_pmids))
    print('Total articles: ', len(pubmed_pmids))

    all_negative_pmids = list(pubmed_pmids - all_positive_pmids)
    all_negative_pmids.sort()
    negative_ratio = 5
    if args.negative_pmids_list is None:
        negative_pmids = set(random.sample(all_negative_pmids, int(len(positive_pmids) * negative_ratio)))
    else:
        print("Loading negatives from list")
        with open(args.negative_pmids_list) as f:
            negative_pmids = set(json.load(f))

    gwas_pmids = gwas.get_gwas_pmids() & pubmed_pmids
    print('Num gwas', len(gwas_pmids))
    clinvar_pmids = clinvar.Clinvar().get_labeled_pmids() & pubmed_pmids
    print('Num clinvar', len(clinvar_pmids))
    all_pmids = list(positive_pmids | negative_pmids | gwas_pmids | clinvar_pmids)
    print('Total papers to download', len(all_pmids))

    random.shuffle(all_pmids)

    with open(args.out_dir + '/dataset_meta.json', 'w') as outfile:
        json.dump({
            'downloaded_pmids': list(all_pmids),
            'positive_pmids': list(positive_pmids),
            'negative_pmids': list(negative_pmids),
            'gwas_pmids': list(gwas_pmids),
            'clinvar_pmids': list(clinvar_pmids),
        }, outfile)

    sys.exit(0)

    articles_per_job = 30
    parts = [all_pmids[i:i+articles_per_job] for i in range(0, len(all_pmids), articles_per_job)]

    input_texts = [
        args.paper_dir + ' ' +
        ' '.join(str(a) for a in pmids) + ' ' +
        '--elsevier_key ' + args.elsevier_key + ' ' + 
        '--sfx_server ' + args.sfx_server
        for pmids in parts]

    # intentionally not localized currently
    results = parasol.submit_jobs(args.conda_path, 'download_papers', input_texts, max_jobs=50)

    parsed_results = [json.loads(result) for result in results]

    print('parsed')

    final_results = {}

    for result in parsed_results:
        final_results.update(result)

    with open(args.out_dir + '/download_meta.json', 'w') as outfile:
        json.dump(final_results, outfile)

    # valid_pmids = {int(pmid) for pmid, pmid_status in final_results.items() if pmid_status['is_valid']}
        
    for pmid, status in final_results.items():
        print(pmid, ': ', status)
