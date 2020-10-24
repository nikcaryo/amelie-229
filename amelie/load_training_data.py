import argparse
import json
import multiprocessing
import os
import pickle
import random
import sys
from collections import defaultdict
import importlib
import constdb

import text_classification

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


def convert_all_to_text(processed_dir, pmids, replace_phenos_with_nothing=False):
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

def shuffle_articles_labels(articles, labels):
    zipped_articles_labels = [x for x in zip(articles, labels)]
    random.shuffle(zipped_articles_labels)
    return [x for x in zip(*zipped_articles_labels)]

def load_training_data(out_dir, process_dir, field_name, replace_phenos_with_nothing=False):
    with open(out_dir + '/dataset_meta.json') as file:
        dataset_info = json.load(file)

    with open(out_dir + '/omim.json') as file:
        omim = json.load(file)
    positives = set(dataset_info['positive_pmids'])
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
    
    print('Good positives: ', len(good_positives), flush=True)
    articles = []
    labels = []

    for pmid, article in convert_all_to_text(process_dir, good_positives,
                                             replace_phenos_with_nothing=replace_phenos_with_nothing):
        omim_data = omim[str(pmid)]
        field = omim_data[field_name]
        articles.append(article)
        labels.append(field[0])

    for pmid, article in convert_all_to_text(process_dir, dataset_info['gwas_pmids'],
                                             replace_phenos_with_nothing=replace_phenos_with_nothing):
        articles.append(article)
        labels.append('gwas')

    articles, labels = shuffle_articles_labels(articles, labels)
    
    
    articles = articles
    labels = labels

    print('Have ', len(articles), flush=True)
    counters = defaultdict(int)
    for label in labels:
        counters[label] += 1

    print('Counts per label:', counters, flush=True)
    print('Done converting', flush=True)
    return articles, labels
#     classifier = text_classification.create_model(articles, labels, cross_val=cross_val, l1=l1)

#     if save:
#         print("SAVING TEXT FIELD RELEVANCE CLASSIFIER", flush=True)
#         with open(out_dir + '/text_field_{}.pkl'.format(field_name), 'wb') as out_file:
#             pickle.dump(classifier, out_file)
#     else:
#         print("NOT SAVING TEXT FIELD RELEVANCE CLASSIFIER", flush=True)
