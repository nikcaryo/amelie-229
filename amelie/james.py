import sys
import sklearn.feature_extraction.text
import sklearn.linear_model
import sklearn.pipeline
import sklearn.ensemble
import sklearn.feature_selection
import sklearn.naive_bayes
from sklearn.decomposition import TruncatedSVD
from sklearn.model_selection import cross_validate
from sklearn.metrics.scorer import make_scorer
from sklearn.metrics import precision_score, recall_score
from sklearn.metrics import classification_report
from load_training_data import load_training_data

import matplotlib.pyplot as plt
from sklearn.metrics import plot_roc_curve
import pandas as pd
from sklearn.metrics import roc_curve, auc

COUNT = 0

def plot_multiclass_roc(clf, X_test, y_test, n_classes, figsize=(17, 6)):
    y_score = clf.decision_function(X_test)

    # structures
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # calculate dummies once
    y_test_dummies = pd.get_dummies(y_test, drop_first=False).values
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test_dummies[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # roc for each class
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot([0, 1], [0, 1], 'k--')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('Receiver operating characteristic example')
    for i in range(n_classes):
        ax.plot(fpr[i], tpr[i], label='ROC curve (area = %0.2f) for label %i' % (roc_auc[i], i))
    ax.legend(loc="best")
    ax.grid(alpha=.4)
    plt.savefig('multiclass_roc_{}.png'.format(COUNT))



def five_fold_cross_val(pipeline, articles, labels):
    print("Number of articles: {}".format(len(articles)))
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
    # featurizers[feat[0]] = feat[1]
    
    # TODO: more shit
    classifiers = {
        #newer logreg args: solver='lbfgs', multi_class='auto'
        "logreg": sklearn.linear_model.LogisticRegression(solver='liblinear',multi_class='ovr',penalty=penalty, max_iter=max_iter),
        "svm": sklearn.svm.LinearSVC(),
        "rfc": sklearn.ensemble.RandomForestClassifier(n_estimators=10), #could be 10, or in newer versions the default is 100
        "nb": sklearn.naive_bayes.MultinomialNB()
    }
    
    # clf = classifiers[model]
    pipe = [('tfidf', sklearn.feature_extraction.text.TfidfVectorizer(min_df=0.01, max_df=0.95))]
    if feat:
        pipe.append(('featurizer', feat))
    pipe.append(('Classifier', classifiers[model]))
    # featurizer = featurizers[feat]
    # featurizer.append(feat)
    # tsvd = TruncatedSVD(n_components=100)
    # rfe = sklearn.feature_selection.RFE(clf,step=1000,verbose=1)
    pipeline = sklearn.pipeline.Pipeline(pipe)
    # pipeline = sklearn.pipeline.Pipeline([
    #     ('Featurizer', featurizer),
    #     # ('RFE', rfe),
    #     # ('tsvd', tsvd),
    #     ('Classifier', clf)
    # ])

    print('Pipeline created', pipeline)
    if cross_val:
        print('Five-fold cross validation')
        five_fold_cross_val(pipeline, articles, labels)
    print('Fitting classifier')
    pipeline.fit(articles, labels)
    # plot_roc_curve(pipeline, articles, labels)
    # plt.show()
    # plot_multiclass_roc(pipeline, articles, labels, n_classes=len(set(labels)))

    return pipeline

if __name__ == "__main__":
    out_dir = "amelie_out_dir"
    process_dir = "amelie_process_dir"
    
    
    # inheritance_modes or variant_types
    # mode = str(sys.argv[1])
    for mode in ["inheritance_modes", "variant_types"]:
        # COUNT += 1
        
        # real deal
        # articles, labels = load_training_data(out_dir, process_dir, mode, new=True)
        
        # small
        articles, labels = load_training_data(out_dir, process_dir, mode, limit=100, new=True)
        article_train, article_test, label_train, label_test = train_test_split(articles, labels, test_size=0.20, random_state=42)

        classifiers = {
            #newer logreg args: solver='lbfgs', multi_class='auto'
            "svm_l2_sqhinge": sklearn.svm.LinearSVC(penalty='l2', loss='squared_hinge'),
            "svm_l1_sqhinge": sklearn.svm.LinearSVC(penalty='l1', loss='squared_hinge'),
            "svm_l2_hinge": sklearn.svm.LinearSVC(penalty='l2', loss='hinge'),
            "svm_l1_hinge": sklearn.svm.LinearSVC(penalty='l1', loss='hinge'),
            "svm_poly": sklearn.svm.SVC(kernel="poly"),
            "svm_rbf": sklearn.svm.SVC(kernel="rbf"),
            "svm_sigmoid": sklearn.svm.SVC(kernel="sigmoid"),
            
        }

        for clf in classifiers.keys():

            featurizers = [
                None,
                sklearn.feature_selection.RFE(classifiers[clf],n_features_to_select=1000,step=1000,verbose=0),
                sklearn.feature_selection.RFE(classifiers[clf],n_features_to_select=5000,step=1000,verbose=0)
            ]
            for feat in featurizers:
                print("!!!!!!!!!!!!!!!!!!! Creating model for mode {}, clf {}, featurizer tf-idf + {}:".format(mode, clf, feat))
                model = create_model(article_train, label_train, clf, feat, cross_val=True)
                predictions = model.predict(article_test)
                pkl_filename = "results/{}-{}-{}-model.pkl".format(mode, clf, feat)
                pickle.dump(model,open(pkl_filename,"wb"))
                #model = pickle.load(model)
        
        # should save this into a text file
                print("For mode {}, clf {}, featurizer tf-idf + {}:".format(mode, clf, feat))
                # print(classification_report(predictions, labels))
                print(classification_report(label_test, predictions))
                print("=========================================")



