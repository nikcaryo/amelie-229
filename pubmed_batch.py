import argparse
import pubMunch

import urllib.request
import gzip
import io
import json
import ftplib

import xml.etree.ElementTree as ET

from ftplib import FTP
import re

import multiprocessing
import constdb

from socket import timeout

path = 'ftp://ftp.ncbi.nlm.nih.gov/pubmed/{}/pubmed{}n{}.xml.gz'

pubmed_xml_regex = r'pubmed([0-9]{2})n([0-9]+).xml.gz$'


def get_num_pieces_and_year(pubmed_type='baseline'):
    for _ in range(5):
        try:
            ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=60)
            ftp.login()
            ftp.cwd('pubmed/' + pubmed_type)
            filenames = ftp.nlst()
            filenames = [x for x in filenames if re.match(pubmed_xml_regex, x)]
            filenames.sort()

            latest_file = filenames[-1]

            print(latest_file)

            version = re.match(pubmed_xml_regex, latest_file)
            return int(version.group(2)), int(version.group(1))
        except timeout:
            print("Socket connection timed out")
    print("Socket timed out 5 times, exiting")
    return None


def process_pubmed_batch_file(pubmed_batch_id, year, pubmed_type='baseline'):
    tree = download_pubmed_batch_file(pubmed_batch_id, year, pubmed_type=pubmed_type)
    return parse_pubmed_xml_file(tree)


def download_pubmed_batch_file(pubmed_batch_id, year, pubmed_type='baseline'):
    year_str = str(year)
    id_str = str(pubmed_batch_id).zfill(4)

    last_exception = None
    # 5 attempts
    for i in range(5):
        try:
            with urllib.request.urlopen(path.format(pubmed_type, year_str, id_str)) as pubmed_file:
                with gzip.open(pubmed_file) as pubmed_xml:
                    return ET.parse(pubmed_xml)
        except Exception as e:
            last_exception = e
    else:
        if last_exception is not None:
            raise last_exception

    return None


def get_abstract_text(abstract):
    return "".join(
        ''.join(a.itertext()) for a in abstract.findall('AbstractText') if a.itertext()
    )


def convert_month(month):
    try:
        rv = int(month)
        return rv
    except ValueError:
        month = month.lower()[:3]
        month_conversion = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
        return month_conversion.index(month) + 1


def get_pmid_and_abstract(pubmed_article):
    citation = pubmed_article.find('MedlineCitation')
    pmid_entry = citation.find('PMID')
    pmid = ''.join(pmid_entry.itertext())
    article = citation.find('Article')

    journal = article.find('Journal')

    journal_title = journal.find('Title')

    if journal_title is not None:
        journal_title_text = ''.join(journal_title.itertext())
    else:
        journal_title_text = None

    version = pmid_entry.get('Version')

    pubmed_data = pubmed_article.find('PubmedData')
    author_list = article.find("AuthorList")
    authors = []
    if author_list is not None:
        for author in author_list.findall('Author'):
            if author.find("LastName") is not None:
                if author.find("LastName").itertext():
                    author_name = ''.join(author.find("LastName").itertext())
                    if author.find("Initials") is not None:
                        if author.find("Initials").itertext():
                            author_name += " " + ''.join(author.find("Initials").itertext())
                authors.append(author_name)

    publish_date = None

    if pubmed_data is not None:
        pubmed_history = pubmed_data.find('History')

        if pubmed_history is not None:
            for pubdate in pubmed_history.findall('PubMedPubDate'):
                if pubdate.get('PubStatus') == 'pubmed':
                    year, month, day = (''.join(pubdate.find('Year').itertext()),
                                        ''.join(pubdate.find('Month').itertext()),
                                        ''.join(pubdate.find('Day').itertext()))

                    publish_date = '{}-{}-{}'.format(year, month, day)

    if publish_date is None:
        journal = article.find('Journal')
        if journal is not None:
            journal_issue = journal.find('JournalIssue')
            if journal_issue is not None:
                pubdate = journal_issue.find('PubDate')
                if pubdate is not None:
                    year = pubdate.find("Year")
                    if year is not None:
                        year = ''.join(year.itertext())
                    month = pubdate.find("Month")
                    if month is not None:
                        month = ''.join(month.itertext())
                        month = convert_month(month)
                    day = pubdate.find("Day")
                    if day is not None:
                        day = ''.join(day.itertext())

                    if year is not None and month is not None and day is not None:
                        publish_date = '{}-{}-{}'.format(year, month, day)

    if version != '1':
        print('Dropping other version of ', pmid)
        return None

    def get_abstract():
        article_abstract = article.find('Abstract')
        other_abstract = citation.find('OtherAbstract')

        if article_abstract is not None:
            return get_abstract_text(article_abstract)
        elif other_abstract is not None:
            return get_abstract_text(other_abstract)
        else:
            return None

    abstract = get_abstract()
    title = ''.join(article.find('ArticleTitle').itertext())

    return pmid, {
        'pmid': pmid,
        'abstract': abstract,
        'title': title,
        'publish_date': publish_date,
        'journal': journal_title_text,
        'authors': authors,
    }


def parse_pubmed_xml_file(tree):
    result = {}

    for article in tree.findall('PubmedArticle'):
        one_article = get_pmid_and_abstract(article)
        if one_article is not None:
            result[one_article[0]] = one_article[1]
    return result

