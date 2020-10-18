
import argparse

from . import parasol, pheno_extractor, gene_extractor, variant_extractor, local_executor
from .refseq_extractor import RefseqExtractor
from . import position_extractor, text
from .dataloaders import pubtator

import json
import constdb
import pickle
import random
import os
import traceback

import subprocess

import sys
import shutil
import tempfile


def exception_to_str(ex):
    return ''.join(traceback.format_exception(etype=type(ex), value=ex, tb=ex.__traceback__))


def process_part_of_paper(source, source_text, extractors, pmid, input_location=None):
    if source_text is None:
        source_text = ''

    orig_text = source_text
    source_text = text.normalize(source_text)
    
    # print('processing', source, len(source_text))
    if 'g_extractor' in extractors:
        gene_mentions = extractors['g_extractor'].extract(source_text, pmid)
    else:
        gene_mentions = 'NO GENE EXTRACTOR DEFINED'
    if 'p_extractor' in extractors:
        phenotype_mentions = list(extractors['p_extractor'].extract(source_text))
    else:
        phenotype_mentions = 'NO PHENO EXTRACTOR DEFINED'
    
    current_result = {
        'raw_text': source_text,
        'orig_text': orig_text,
        'gene_mentions': gene_mentions,
        'phenotype_mentions': phenotype_mentions
    }

    if source == 'text':
        if 'refseq_extractor' in extractors and 'v_extractor' in extractors and 'pos_extractor' in extractors:
            ensembl_to_mentioned_refseq = extractors['refseq_extractor'].extract(orig_text)
            variant_mentions = list(extractors['v_extractor'].extract(orig_text,
                                                                      gene_mentions,
                                                                      ensembl_to_mentioned_refseq))

            strings_to_extract = set()
            for mention in variant_mentions:
                strings_to_extract.add(mention.variant_text)
                for modification in mention.modifications:
                    strings_to_extract.add(modification.grounded_refseq.split('.')[0])
                    strings_to_extract.add(modification.nm_refseq.split('.')[0])

            all_gene_mentions = {mention for mentions in gene_mentions.values() for mention in mentions}

            for mention in all_gene_mentions:
                strings_to_extract.update(mention.orig_words.split('|^|'))

            if len(strings_to_extract) != 0 and orig_text != '' and input_location is not None and \
                    os.path.isfile(input_location) and os.stat(input_location).st_size > 0:
                positions = extractors['pos_extractor'].extract(input_location, list(strings_to_extract))
            else:
                positions = {}

            current_result.update({
                'variant_mentions': list(variant_mentions),
                'positions': positions
            })
        else:
            print("Cannot extract variants; at least one of refseq_extractor, v_extractor or pos_extractor is missing;"
                  "adding empty variants list")
            current_result.update({
                'variant_mentions': list(),
                'positions': {}
            })

    return current_result


def get_text_for_paper(input_location, pubmed_info, use_ta_if_download_failed=False):
    if os.path.isfile(input_location) and os.stat(input_location).st_size > 0:
        result = subprocess.run(['pdftotext', '-enc', 'UTF-8', input_location, '-'],  stdout=subprocess.PIPE)

        result_text = result.stdout.decode('utf8', 'replace').strip()

        if len(result_text) > 500000:
            # This is probably junk
            result_text = ''
    elif use_ta_if_download_failed:
        result_text = str(pubmed_info['title']) + ' . ' + str(pubmed_info['abstract'])
    else:
        result_text = None
        assert False, "Don't call this method unless you have some text"

    return result_text


def process_paper(paper_dir, process_dir, extractors, pubmed_info, pmid,
                  fulltext_only=False, use_ta_if_download_failed=False):
    input_location = paper_dir + '/' + str(pmid) + '.pdf'
    output_location = process_dir + '/' + str(pmid) + '.pkl'

    if os.path.isfile(output_location):
        print("Already have output %s, returning valid" % output_location)
        return {
            'is_valid': True,
            'status': 'Already processed'
        }

    print("Input location: {}".format(input_location))
    if not os.path.isfile(input_location) or os.stat(input_location).st_size == 0:
        if not use_ta_if_download_failed:
            print("Do not have input %s or is empty" % input_location)
            return {
                'is_valid': False,
                'status': 'Not downloaded'
            }

    try:
        if fulltext_only:
            text_sources = {
                'text': get_text_for_paper(input_location,
                                           pubmed_info=pubmed_info,
                                           use_ta_if_download_failed=False),
            }
        else:
            text_sources = {
                'abstract': pubmed_info['abstract'],
                'title': pubmed_info['title'],
                'text': get_text_for_paper(input_location,
                                           pubmed_info=pubmed_info,
                                           use_ta_if_download_failed=use_ta_if_download_failed),
            }

        text_results = {}

        for source, source_text in text_sources.items():
            text_results[source] = process_part_of_paper(source, source_text, extractors, pmid, input_location)
        # text_results['meta'] = {"publish_date": pubmed_info['publish_date'], "journal": pubmed_info['journal']}

        print("Dumping output %s" % output_location)
        with open(output_location, 'wb') as output:
            pickle.dump(text_results, output)

        return {
            'is_valid': True,
            'status': 'Processed'
        }

    except Exception as e:
        exception_string = exception_to_str(e)
        print("Have exception %s, returning error" % exception_string)
        return {
            'is_valid': False,
            'status': 'Error',
            'error_message': exception_to_str(e)
        }


def create_pdf(dir, name, text): 
    input_text_file = tempfile.NamedTemporaryFile(mode='w')
    print(text, file=input_text_file)
    input_text_file.flush()
    enscript = subprocess.Popen(['/usr/bin/enscript', '--word-wrap',
                                 '-B', input_text_file.name, '-o', '-'], stdout=subprocess.PIPE)
    ps2pdf = subprocess.Popen(['/usr/bin/ps2pdf', '-', os.path.join(dir, name)], stdin=enscript.stdout)
    ps2pdf.wait()


def process_title_abstract_only(paper_dir, process_dir, extractors, pubmed_info, pmid):
    output_location = process_dir + '/' + str(pmid) + '.pkl'

    if os.path.isfile(output_location):
        print("Already have output %s, returning valid" % output_location)
        return {
            'is_valid': True,
            'status': 'Already processed'
        }

    # next lines only useful if using variant extractor, but the variant extractor doesn't work anyways ...
    # pdf_name = '%s.pdf' % pmid
    # if not os.path.isfile(os.path.join(paper_dir, pdf_name)):
    #     create_pdf(dir=paper_dir, name=pdf_name, text=str(pubmed_info['title'])
    # + '\n\n' + str(pubmed_info['abstract']))

    try:
        text_sources = {
            'abstract': pubmed_info['abstract'],
            'title': pubmed_info['title'],
            'text': str(pubmed_info['title']) + ' ' + str(pubmed_info['abstract'])
        }

        text_results = {}

        for source, source_text in text_sources.items():
            text_results[source] = process_part_of_paper(source, source_text, extractors, pmid)
        # text_results['meta'] = {"publish_date": pubmed_info['publish_date'], "journal": pubmed_info['journal']}

        print("Dumping output %s" % output_location)
        with open(output_location, 'wb') as output:
            pickle.dump(text_results, output)

        print("Returning valid")
        return {
            'is_valid': True,
            'status': 'Processed'
        }

    except Exception as e:
        exception_string = exception_to_str(e)
        print("Have exception %s, returning error" % exception_string)
        return {
            'is_valid': False,
            'status': 'Error',
            'error_message': exception_string
        }


def process_papers_program():
    parser = argparse.ArgumentParser()

    parser.add_argument('out_file', 
        type=str,
        help='The name of the database to output')

    parser.add_argument('out_dir', 
        type=str,
        help='The name of the database to output')

    parser.add_argument('paper_dir',
                        type=str)

    parser.add_argument('process_dir', 
        type=str, 
        default='processed_papers',
        help='The name of the database to output')

    parser.add_argument('pubmed_ids', 
        type=int,
        nargs='*',
        help='The pubmed ids to download')

    parser.add_argument('--pubmunch_data_dir', 
        type=str, 
        default='/cluster/u/ethanid/pubmunchData', 
        help='The location of your pubmed baseline dump'
    )

    parser.add_argument('--title_abstract_only',
                        action="store_true",
                        default=False)

    args = parser.parse_args()
    print(args, file=sys.stderr)

    if not os.path.isdir(args.out_dir):
        print("Out dir %s does not exist" % args.out_dir)
        sys.exit(1)

    if not args.title_abstract_only:
        if not os.path.isdir(args.paper_dir):
            print("Paper dir %s does not exist" % args.paper_dir)
            sys.exit(1)

    if not os.path.isdir(args.pubmunch_data_dir):
        print("Pubmunch data dir %s does not exist" % args.pubmunch_data_dir)
        sys.exit(1)

    with pubtator.Pubtator(args.out_dir + '/pubtator.db') as pub:
        with position_extractor.PositionExtractor() as pos_extractor:
            p_extractor = pheno_extractor.PhenoExtractor()
            g_extractor = gene_extractor.GeneExtractor(pub)
            v_extractor = variant_extractor.VariantExtractor(args.pubmunch_data_dir)
            refseq_extractor = RefseqExtractor()

            extractors = {
                'p_extractor': p_extractor,
                'g_extractor': g_extractor,
                'v_extractor': v_extractor,
                'pos_extractor': pos_extractor,
                'refseq_extractor': refseq_extractor,
            }

            print('About to read stored pubmed', file=sys.stderr)
            with constdb.read(args.out_dir + '/pubmed.db', mmap=False, keys_to_read=set(args.pubmed_ids)) as pubmed_dump:
                print('Ready ', args.pubmed_ids, file=sys.stderr)
                results = {}

                for pmid in args.pubmed_ids:
                    dumped = pubmed_dump.get(pmid)
                    if dumped is None:
                        continue
                    pubmed_info = json.loads(dumped.decode('utf8'))
                    if args.title_abstract_only:
                        results[pmid] = process_title_abstract_only(
                            args.paper_dir,
                            args.process_dir,
                            extractors,
                            pubmed_info,
                            pmid
                        )
                    else:
                        results[pmid] = process_paper(
                            args.paper_dir,
                            args.process_dir,
                            extractors,
                            pubmed_info,
                            pmid
                        )

    print(results)

    with open(args.out_file, 'w') as outfile:
        json.dump(results, outfile)


def process_all_papers_program():
    parser = argparse.ArgumentParser()

    parser.add_argument('out_dir', 
        type=str,
        help='The name of the database to output')

    parser.add_argument('paper_dir',
        type=str, 
        help='The name of the database to output')

    parser.add_argument('process_dir', 
        type=str, 
        default='processed_papers',
        help='The name of the database to output')

    parser.add_argument('--conda_path', 
        type=str, 
        default='/cluster/u/jbirgmei/miniconda3/bin')

    parser.add_argument('--pubmunch_data_dir', 
        type=str, 
        default='/cluster/u/jbirgmei/shallowfloat/amelie_ethan/pubmunchData')

    parser.add_argument('--title_abstract_only',
                        action="store_true",
                        default=False)

    parser.add_argument('--process_only_first',
                        type=int,
                        default=None)

    parser.add_argument('--delete_parasol_dir',
                        action="store_true",
                        default=False)

    parser.add_argument('--local_processing',
                        action='store_true',
                        default=False)

    args = parser.parse_args()
    print(args)

    if not os.path.isdir(args.out_dir):
        print("Out dir %s does not exist" % args.out_dir)
        sys.exit(1)

    if not os.path.isdir(args.process_dir):
        print("Process dir %s does not exist" % args.process_dir)
        sys.exit(1)

    if not args.title_abstract_only:
        if not os.path.isdir(args.paper_dir):
            print("Paper dir %s does not exist" % args.paper_dir)
            sys.exit(1)

    if not os.path.isdir(args.pubmunch_data_dir):
        print("Pubmunch data dir %s does not exist" % args.pubmunch_data_dir)
        sys.exit(1)

    # os.mkdir(args.process_dir)  # fool. This cost me 20 hrs

    with open(args.out_dir + '/dataset_meta.json') as input_meta_file:
        input_meta = json.load(input_meta_file)

    print(input_meta.keys())

    pmids_to_process = list(set(input_meta['clinvar_pmids'] + input_meta['gwas_pmids'] +
                                input_meta['negative_pmids'] + input_meta['positive_pmids']))
    pmids_to_process.sort()

    random.seed(45123678451)
    random.shuffle(pmids_to_process)

    if args.process_only_first is not None:
        pmids_to_process = pmids_to_process[:args.process_only_first]
    print("All PMIDs: %d" % len(pmids_to_process))

    articles_per_job = 50

    not_processed_pmids = []
    for pmid in pmids_to_process:
        output_location = args.process_dir + '/' + str(pmid) + '.pkl'
        if not os.path.isfile(output_location):
            not_processed_pmids.append(pmid)
    pmids_to_process = not_processed_pmids
    print("Not processed PMIDs: %d" % len(pmids_to_process))

    if not args.title_abstract_only:
        downloaded_pmids_to_process = []
        for pmid in pmids_to_process:
            input_location = args.paper_dir + '/' + str(pmid) + '.pdf'
            if os.path.isfile(input_location) and os.stat(input_location).st_size > 0:
                downloaded_pmids_to_process.append(pmid)
        pmids_to_process = downloaded_pmids_to_process
        print("Downloaded PMIDs to process: %d" % len(pmids_to_process))

    print('Going to process', len(pmids_to_process))
    parts = [pmids_to_process[i:i+articles_per_job] for i in range(0, len(pmids_to_process), articles_per_job)]

    if args.title_abstract_only:
        input_texts = [
            '--pubmunch_data_dir ' + args.pubmunch_data_dir + ' ' +
            '--title_abstract_only ' +
            args.out_dir + ' ' +
            args.paper_dir + ' ' +
            args.process_dir + ' ' +
            ' '.join(str(a) for a in pmids)
            for pmids in parts
        ]
    else:
        input_texts = [
            '--pubmunch_data_dir ' + args.pubmunch_data_dir + ' ' +
            args.out_dir + ' ' +
            args.paper_dir + ' ' +
            args.process_dir + ' ' +
            ' '.join(str(a) for a in pmids)
            for pmids in parts
        ]

    print(input_texts[:3])

    if args.local_processing:
        parasol_result = local_executor.submit_jobs(args.conda_path, 'process_papers',
                                                    input_texts)
    else:
        parasol_result = parasol.submit_jobs(args.conda_path,
                                             'process_papers',
                                             input_texts, max_jobs=200)
    print(parasol_result['output'])
    parsed_results = [json.loads(output) for output in parasol_result['output']]
    print('parsed')

    final_results = {}

    for result in parsed_results:
        final_results.update(result)

    with open(args.out_dir + '/process_meta.json', 'w') as outfile:
        json.dump(final_results, outfile)
        
    for pmid, status in final_results.items():
        print(pmid, ': ', status)

    if args.delete_parasol_dir:
        shutil.rmtree(parasol_result['run_directory'])

