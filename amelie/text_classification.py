
import nltk
from nltk.corpus import stopwords
import sklearn.feature_extraction.text
import sklearn.linear_model
import sklearn.pipeline
from sklearn.model_selection import cross_validate
from sklearn.metrics.scorer import make_scorer
from sklearn.metrics import precision_score, recall_score

from dataloaders.hpo import HPO

stopwords = set(stopwords.words('english'))
print(__name__)
hpo = HPO()


def flatten(l):
    return [item for sublist in l for item in sublist]


def remove_stopwords(words):
    return [w for w in words if not w in stopwords]


def replace_genes(words, gms):
    for gm in gms:
        for i in range(len(gm.wordidxs_start)):
            start = gm.wordidxs_start[i]
            end = gm.wordidxs_end[i]
            if start > len(words):
                print(start, end)
                print(words)
            # assert start < len(words), (start, gm, words, len(words), pmid)
            words[start] = 'XGENE'
            for i in range(start + 1, end):
                words[i] = ''
    return words


def replace_phenos(words, pms, bad_subtrees, replace_phenos_with_nothing):
    for pm in pms:
        replace = True
        for bad_subtree in bad_subtrees:
            if hpo.has_descendant(bad_subtree, pm.hpo_id):
                replace = False
                break
        if not replace:
            continue
        wordidxs = pm.wordidxs
        # assert wordidxs[0] < len(words), (wordidxs[0], pm, words, len(words), pmid)
        if replace_phenos_with_nothing:
            words[wordidxs[0]] = ''
        else:
            words[wordidxs[0]] = 'XPHENO'
        for i in wordidxs[1:]:
            words[i] = ''
    return words


def process_text(processed_data, prefix, replace_phenos_with_nothing):
    all_gene_mentions = {mention for mentions in processed_data['gene_mentions'].values() for mention in mentions}
    words = [word.lower() for word in processed_data['raw_text'].split()]
    words = replace_genes(words, all_gene_mentions)
    bad_subtrees = ['HP:0002664', 'HP:0012125', 'HP:0100753', 'HP:0003254']
    words = replace_phenos(words, processed_data['phenotype_mentions'], bad_subtrees=bad_subtrees,
                           replace_phenos_with_nothing=replace_phenos_with_nothing)
    words = remove_stopwords(words)
    if replace_phenos_with_nothing:
        return ' '.join(prefix + word for word in words if len(word) > 0)
    else:
        return ' '.join(prefix + word for word in words)


def convert_to_text(processed_article, replace_phenos_with_nothing, use_main_text):
    result_data = process_text(processed_article['title'], 'TITLE_',
                               replace_phenos_with_nothing=replace_phenos_with_nothing)
    if processed_article['abstract'] is not None:
        result_data += ' ' + process_text(processed_article['abstract'], 'ABSTRACT_',
                                          replace_phenos_with_nothing=replace_phenos_with_nothing)
    if use_main_text:
        result_data += ' ' + process_text(processed_article['text'], 'TEXT_',
                                          replace_phenos_with_nothing=replace_phenos_with_nothing)
    return result_data


def five_fold_cross_val(pipeline, articles, labels):
    num_labels = len(set(labels))

    if num_labels <= 2:
        print("Scoring: binary")
        scoring = {'precision': make_scorer(precision_score, average='binary'),
                   'recall': make_scorer(recall_score, average='binary')}
    else:
        print("Scoring: micro")
        scoring = {'precision': make_scorer(precision_score, average='micro'),
                   'recall': make_scorer(recall_score, average='micro')}

    scores = sklearn.model_selection.cross_validate(pipeline,
                                                    articles,
                                                    labels,
                                                    cv=5,
                                                    verbose=1,
                                                    scoring=scoring)

    print("Scores: ")
    print(scores)


def create_model(articles, labels, cross_val=False, l1=False):
    if l1:
        clf = sklearn.linear_model.LogisticRegression(penalty='l1', max_iter=1000)
    else:
        clf = sklearn.linear_model.LogisticRegression(penalty='l2', max_iter=1000)
    pipeline = sklearn.pipeline.Pipeline([
        ('Featurizer', sklearn.feature_extraction.text.TfidfVectorizer(min_df=0.01, max_df=0.95)),
        ('Classifier', clf)
    ])

    print('Pipeline created')
    if cross_val:
        print('Five-fold cross validation')
        five_fold_cross_val(pipeline, articles, labels)
    print('Fitting classifier')
    pipeline.fit(articles, labels)

    return pipeline
