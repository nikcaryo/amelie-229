import sklearn.feature_extraction.text
import sklearn.linear_model
import sklearn.pipeline
from sklearn.model_selection import cross_validate
from sklearn.metrics.scorer import make_scorer
from sklearn.metrics import precision_score, recall_score

from load_training_data import load_training_data


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

if __name__ == "__main__":
    out_dir = "amelie_out_dir"
    process_dir = "amelie_process_dir"

    articles, labels = load_training_data(out_dir, process_dir, "inheritance_modes", num=10)
    print(labels)
    
    create_model(articles, labels, cross_val=True)





