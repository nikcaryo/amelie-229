
import argparse
import sys, os
import sqlite3
import glob
from collections import defaultdict


def ta_stats(directory):
    num_ta = 0
    relevant_pmid_titles_journals = {}
    for filename in glob.glob(directory + "/out*"):
        with open(filename) as f:
            for line in f:
                if line.startswith('Tested title/abstract of PMID'):
                    num_ta += 1
                    continue
                if not line.startswith('After testing for relevance'):
                    continue
                line = line.strip().split('\t')
                # have stupid stdout clashes in rare cases, this is a way of recovering like 99% of the good stuff:
                pmid, title, journal = line[1], line[2], line[3]
                relevant_pmid_titles_journals[pmid] = {'title': title, 'journal': journal}
    return {'num_ta': num_ta, 'relevant_pmid_titles_journals': relevant_pmid_titles_journals}


def get_downloaded(paper_dir, pmids):
    rv = set()
    for pmid in pmids:
        input_location = os.path.join(paper_dir, pmid + '.pdf')
        if os.path.isfile(input_location) and os.stat(input_location).st_size > 0:
            rv.add(pmid)
    return rv


def get_loaded(amelie_knowledgebase, pmids):
    connection = sqlite3.connect(amelie_knowledgebase)
    cur = connection.cursor()
    loaded_pmids = set(x[0] for x in cur.execute('select distinct pmid from amelie_paper').fetchall())
    return loaded_pmids


def get_top_stats(relevant_pmid_titles_journals, downloaded_pmids):
    num_journal_downloads = defaultdict(lambda: 0)
    num_journal_fails = defaultdict(lambda: 0)
    for pmid in relevant_pmid_titles_journals:
        journal = relevant_pmid_titles_journals[pmid]['journal']
        if pmid in downloaded_pmids:
            num_journal_downloads[journal] += 1
        else:
            num_journal_fails[journal] += 1
    top_downloaded = sorted((y, x) for x, y in num_journal_downloads.items())[::-1]
    top_failed = sorted((y, x) for x, y in num_journal_fails.items())[::-1]
    return top_downloaded, top_failed


def extraction_report():
    parser = argparse.ArgumentParser()
    parser.add_argument('extractor_tmpdir', type=str)
    parser.add_argument('amelie_paper_dir', type=str)
    parser.add_argument('amelie_knowledgebase', type=str)
    args = parser.parse_args()

    rv = ta_stats(args.extractor_tmpdir)
    num_ta = rv['num_ta']
    relevant_pmid_titles_journals = rv['relevant_pmid_titles_journals']
    downloaded_pmids = get_downloaded(paper_dir=args.amelie_paper_dir, pmids=relevant_pmid_titles_journals)
    loaded_pmids = get_loaded(amelie_knowledgebase=args.amelie_knowledgebase, pmids=downloaded_pmids)

    print("Processed %d total new titles and abstracts" % num_ta)
    print("Possibly relevant: %d" % len(relevant_pmid_titles_journals))
    print("Successfully downloaded: %d" % len(downloaded_pmids))
    print("Successfully loaded into database: %d" % len(loaded_pmids.intersection(downloaded_pmids)))
    print()

    top_downloaded, top_failed = get_top_stats(relevant_pmid_titles_journals, downloaded_pmids)
    print("Top downloaded journals:")
    for num, journal in top_downloaded:
        print("\t%d\t%s" % (num, journal))
    print()
    print("Top failed journals:")
    for num, journal in top_failed:
        print("\t%d\t%s" % (num, journal))
    print()

    print()
    print("Relevant articles info:")
    for pmid in relevant_pmid_titles_journals:
        print('Downloaded: %d, Loaded into database: %d\n\tPMID %s, Journal: "%s", Title: "%s"' %
              (pmid in downloaded_pmids,
               pmid in loaded_pmids,
               pmid,
               relevant_pmid_titles_journals[pmid]['journal'],
               relevant_pmid_titles_journals[pmid]['title'],
               ))
