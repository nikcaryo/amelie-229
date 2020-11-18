import sys
import sklearn.feature_extraction.text
import sklearn.linear_model
import sklearn.pipeline
from sklearn import tree
from sklearn.model_selection import cross_validate
from sklearn.metrics.scorer import make_scorer
from sklearn.metrics import precision_score, recall_score
from sklearn.metrics import classification_report
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

def create_model(articles, labels, model, feat, 
                 cross_val=False, 
                 penalty="l2", 
                 max_iter=1000):
    
    # TODO: more shit
    featurizers = {
        "tfidf": sklearn.feature_extraction.text.TfidfVectorizer(min_df=0.01, max_df=0.95)
         
    }
    
    # TODO: more shit
    classifiers = {
        "logreg": sklearn.linear_model.LogisticRegression(penalty=penalty, max_iter=max_iter),
        "svm": sklearn.svm.LinearSVC(),
        "tree": tree.DecisionTreeClassifier()
    }
    
    clf = classifiers[model]
    featurizer = featurizers[feat]
    
    pipeline = sklearn.pipeline.Pipeline([
        ('Featurizer', featurizer),
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
    
    
    # inheritance_modes or variant_types
    # mode = str(sys.argv[1])
    for mode in ["inheritance_modes", "variant_types"]:
        
        articles, labels = load_training_data(out_dir, process_dir, mode, new=False)
#         articles, labels = load_training_data(out_dir, process_dir, mode, new=True)
#         for clf in ["logreg"]:
#             print(f"MODE: {mode}")
#             print(f"CLF: {clf}")
      
#             model = create_model(articles, labels, "tree", "tfidf", cross_val=True)
#             predictions = model.predict(articles)
        
#             # should save this into a text file
#             print(classification_report(labels, predictions))
#             print("=========================")



