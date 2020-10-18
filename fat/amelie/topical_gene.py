import sklearn.pipeline
import sklearn.feature_extraction.text

from .text_classification import five_fold_cross_val

import numpy as np


class ItemSelector(sklearn.base.TransformerMixin):
    def __init__(self, key):
        self.key = key

    def fit(self, x, y=None):
        return self

    def transform(self, X):
        return [data_dict[self.key] for data_dict in X]


class CustomNumbersExtractor(sklearn.base.TransformerMixin):
    def fit(self, x, y=None):
        return self

    def transform(self, X):
        return [np.array(data) for data in X]


def get_word_windows(prefix, genes):
    rv = []
    for gm in genes:
        for words in gm.surrounding_words_left:
            rv.extend(prefix + '_LEFT_' + x for x in words.split('|^|') if x != "")
        for words in gm.surrounding_words_right:
            rv.extend(prefix + '_RIGHT_' + x for x in words.split('|^|') if x != "")
    return rv


def count_in(gms):
    result = 0
    for gm in gms:
        result += gm.num_mentions
    return result


class TopicalGeneFeaturizer(sklearn.base.TransformerMixin):
    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return [self.transform_article(processed_article, eid) for processed_article, eid in X]

    def transform_article(self, processed_article, eid):
        word_windows = []

        for source, source_items in processed_article.items():
            word_windows.extend(get_word_windows(source, source_items['gene_mentions'][eid]))

        possible_sources = ['title', 'abstract', 'text']
        count_features = []

        for possible_source in possible_sources:
            count_features.append(count_in(processed_article[possible_source]['gene_mentions'][eid]))

        return {
            'words': ' '.join(word_windows),
            'custom': count_features,
        }


def create_model(articles, labels, cross_val=False, l1=False):
    vect = sklearn.feature_extraction.text.TfidfVectorizer(min_df=0.01, max_df=0.95)
    feature_union = sklearn.pipeline.FeatureUnion(transformer_list=
           [('words',
             sklearn.pipeline.Pipeline([('selector',
                                         ItemSelector(key='words')),
                                        ('tfidf', vect)])),
            ('custom',
             sklearn.pipeline.Pipeline([('selector',
                                         ItemSelector(key='custom')),
                                        ('customExtractor',
                                         CustomNumbersExtractor())]))])

    if l1:
        clf = sklearn.linear_model.LogisticRegression(penalty='l1', max_iter=300)
    else:
        clf = sklearn.linear_model.LogisticRegression(penalty='l2', max_iter=300)

    pipeline = sklearn.pipeline.Pipeline([
        ('Topical gene converter', TopicalGeneFeaturizer()),
        ('Featurizer', feature_union),
        ('Classifier', clf)
    ])

    print('Pipeline created')
    if cross_val:
        print('Five-fold cross-validation')
        five_fold_cross_val(pipeline, articles, labels)
    print('Fitting classifier')
    pipeline.fit(articles, labels)

    return pipeline
