
from collections import namedtuple, defaultdict

from math import sqrt, acos

import sklearn.ensemble
import sklearn.pipeline
import sklearn.model_selection

from .text_classification import five_fold_cross_val


def vector_length(v):
    return sqrt(v[0] * v[0] + v[1] * v[1])


def euclidean(p1, p2):
    v = (p2['x'] - p1['x'], p2['y'] - p1['y'])
    return vector_length(v)


def inv_cos_angle(p1, p2):
    v = (p2['x'] - p1['x'], p2['y'] - p1['y'])
    vlen = vector_length(v)
    if vlen == 0:
        return 0
    cosphi = (p2['y'] - p1['y']) / vector_length(v)
    return acos(cosphi)


RelativePosition2D = namedtuple("RelativePosition2D", ["page_dist", "dist", "angle", "relx",
                                                       "rely", "above", "lstop1", "rstop1",
                                                       "lletter1", "rletter1", "lstop2", "rstop2",
                                                       "lletter2", "rletter2"])
noneRelPos2D = RelativePosition2D(page_dist=0, dist=0, angle=0, relx=0, rely=0, above=float(False),
                                  lstop1=0, rstop1=0, lletter1=0,
                                  rletter1=0, lstop2=0, rstop2=0, lletter2=0, rletter2=0)

RelativePositionWord = namedtuple("RelativePositionWord", ["dist", "above", "lstop1", "rstop1",
                                                           "lletter1", "rletter1", "lstop2", "rstop2",
                                                           "lletter2", "rletter2"])
noneRelPosWord = RelativePositionWord(dist=0, above=float(False),
                                      lstop1=0, rstop1=0, lletter1=0, rletter1=0,
                                      lstop2=0, rstop2=0, lletter2=0, rletter2=0)


def relationship_2d(gpos, vpos):
    page_dist = gpos['page'] - vpos['page']
    dist = euclidean(gpos, vpos)
    angle = inv_cos_angle(gpos, vpos)
    relx = gpos['x'] - vpos['x']
    rely = gpos['y'] - vpos['y']
    if rely <= 20:
        above = float(True)
    else:
        above = float(False)
    # topical gene score of gene is added when function is called
    return RelativePosition2D(page_dist=page_dist, dist=dist, angle=angle,
                              relx=relx, rely=rely, above=above,
                              lstop1=gpos['numStopwordsLeft'], rstop1=gpos['numStopwordsRight'],
                              lletter1=gpos['numLetterLeft'], rletter1=gpos['numLetterRight'],
                              lstop2=vpos['numStopwordsLeft'], rstop2=vpos['numStopwordsRight'],
                              lletter2=vpos['numLetterLeft'], rletter2=vpos['numLetterRight'])


def relationship_wordspace(gpos, vpos):
    # topical gene score of gene is added when function is called
    return RelativePositionWord(dist=abs(gpos['wordIndex'] - vpos['wordIndex']),
                                above=float(gpos['wordIndex'] - vpos['wordIndex'] < 0),
                                lstop1=gpos['numStopwordsLeft'], rstop1=gpos['numStopwordsRight'],
                                lletter1=gpos['numLetterLeft'], rletter1=gpos['numLetterRight'],
                                lstop2=vpos['numStopwordsLeft'], rstop2=vpos['numStopwordsRight'],
                                lletter2=vpos['numLetterLeft'], rletter2=vpos['numLetterRight'])


PdfPosition = namedtuple('PdfPosition', ['pmid', 'name', 'page', 'x', 'y', 'orientation',
                                         'wordidx', 'numStopwordsLeft', 'numStopwordsRight',
                                         'numLetterLeft', 'numLetterRight'])


def get_pdf_positions(pmid, orig_str, cur):
    rv = []
    cur.execute('SELECT * FROM positions WHERE pmid=? and name=?', [pmid, orig_str])
    rows = cur.fetchall()
    for row in rows:
        # print(row)
        # originally named, num_non_letter, the number in reality is num_letter
        pmid, name, page, x, y, orientation, wordidx, num_stopwords_left, \
            num_stopwords_right, num_letter_left, num_letter_right = row
        rv.append(PdfPosition(pmid=pmid, name=name, page=page, x=x, y=y,
                              orientation=orientation,
                              wordidx=wordidx, num_stopwords_left=num_stopwords_left,
                              num_stopwords_right=num_stopwords_right,
                              num_letter_left=num_letter_left,
                              num_letter_right=num_letter_right))
    return rv


def flatten(l):
    return [x for y in l for x in y]


def double_flatten(l):
    return [x for y in l for x in y[0]], [x for y in l for x in y[1]]


def find_closest_above_2d(rel_positions):
    rv_list = [x for x in sorted(rel_positions, key=lambda relpos: relpos.dist) if x.above and x.page_dist == 0]
    feature_names = ["found_closest_above_2d"] + \
                    ["%s_above_2d" % name for name, _ in noneRelPos2D._asdict().items()]
    if rv_list:
        return feature_names, [float(True)] + [x for x in rv_list[0]], rel_positions[rv_list[0]]
    return feature_names, [float(False)] + [x for x in noneRelPos2D], None


def find_closest_above_wordspace(rel_positions):
    rv_list = [x for x in sorted(rel_positions, key=lambda relpos: relpos.dist) if x.above]
    feature_names = ["found_closest_above_wordspace"] + \
                    ["%s_above_wordspace" % name for name, _ in noneRelPosWord._asdict().items()]
    if rv_list:
        return feature_names, [float(True)] + [x for x in rv_list[0]], rel_positions[rv_list[0]]
    return feature_names, [float(False)] + [x for x in noneRelPosWord], None


def find_closest_2d(rel_positions):
    rv_list = [x for x in sorted(rel_positions, key=lambda relpos: relpos.dist) if x.page_dist == 0]
    feature_names = ["found_closest_2d"] + \
                    ["%s_closest_2d" % name for name, _ in noneRelPos2D._asdict().items()]
    if rv_list:
        return feature_names, [float(True)] + [x for x in rv_list[0]], rel_positions[rv_list[0]]
    return feature_names, [float(False)] + [x for x in noneRelPos2D], None


def find_closest_wordspace(rel_positions):
    rv_list = [x for x in sorted(rel_positions, key=lambda relpos: relpos.dist)]
    feature_names = ["found_closest_wordspace"] + \
                    ["%s_closest_wordspace" % name for name, _ in noneRelPosWord._asdict().items()]
    if rv_list:
        return feature_names, [float(True)] + [x for x in rv_list[0]], rel_positions[rv_list[0]]
    return feature_names, [float(False)] + [x for x in noneRelPosWord], None


class DocumentVariantFeatureExtractor:
    def __init__(self, processed_article, topical_gene_classifier):
        self.gene_to_positions = {}
        self.gene_to_topical_score = {}
        self.positions_to_genes = defaultdict(set)

        all_gene_positions = {}

        for gene, mentions in processed_article['text']['gene_mentions'].items():
            probas = topical_gene_classifier.predict_proba([(processed_article, gene)])
            assert probas.shape[0] == 1, probas.shape
            score = probas[0, 1]
            self.gene_to_topical_score[gene] = score

            positions = {}
            for mention in mentions:
                for gene_text in mention.orig_words.split('|^|'):
                    if gene_text in processed_article['text']['positions']:
                        for position in processed_article['text']['positions'][gene_text]:
                            positions[position['wordIndex']] = position
                            all_gene_positions[position['wordIndex']] = position
                            self.positions_to_genes[position['wordIndex']].add(gene)

            self.gene_to_positions[gene] = list(positions.values())

        self.position_map = processed_article['text']['positions']
        self.all_gene_positions = list(all_gene_positions.values())

    def featurize_position(self, variant_mention, vpos):
        result = {}

        distance_funcs = {
            'closest_2d': (relationship_2d, find_closest_2d),
            'closest_above_2d': (relationship_2d, find_closest_above_2d),
            'closest_wordspace': (relationship_wordspace, find_closest_wordspace),
            'closest_above_wordspace': (relationship_wordspace, find_closest_above_wordspace),
        }

        distance_func_to_pos = {}

        for name, (relationship_func, distance_func) in distance_funcs.items():
            relpos_to_gpos_local = {relationship_func(gpos, vpos): gpos
                                    for gpos in self.gene_to_positions[variant_mention.gene]}
            feature_names, local_closest, local_gpos = distance_func(relpos_to_gpos_local)

            relpos_to_gpos_all = {relationship_func(gpos, vpos): gpos for gpos in self.all_gene_positions}
            _, all_closest, all_gpos = distance_func(relpos_to_gpos_all)

            if all_gpos is None:
                right_gene_score = 0
            else:
                all_eids = self.positions_to_genes[all_gpos['wordIndex']]
                right_gene_score = max(self.gene_to_topical_score[eid] for eid in all_eids)

            names = [x + '_target' for x in feature_names] + \
                    [x + '_all' for x in feature_names] + \
                    [name + '_right_gene_score']
            features = local_closest + all_closest + [right_gene_score]
            distance_func_to_pos[name] = local_gpos, all_gpos

            for name, feature in zip(names, features):
                result[name] = feature

        groups = [
            ['closest_2d', 'closest_above_2d'],
            ['closest_wordspace', 'closest_above_wordspace'],
        ]

        for group in groups:
            for item in group:
                for other_item in group:
                    for i, position in enumerate(['local', 'all']):
                        if item == other_item and position == 'local':
                            continue  # No need to compare local to itself

                        feature_name = '{}_local={}_{}'.format(item, other_item, position)
                        feature_value = distance_func_to_pos[item][0] == distance_func_to_pos[other_item][i]
                        result[feature_name] = feature_value

        gene_score = self.gene_to_topical_score[variant_mention.gene]
        result["eid_right_gene_score"] = gene_score
        return result

    def featurize(self, variant_mention):
        results = []
        for vpos in self.position_map.get(variant_mention.variant_text, []):
            results.append(self.featurize_position(variant_mention, vpos))
        return results


def create_model(articles, labels, cross_val=False):
    pipeline = sklearn.pipeline.Pipeline([
        ('Featurizer', sklearn.feature_extraction.DictVectorizer()),
        ('Classifier', sklearn.ensemble.GradientBoostingClassifier(n_estimators=100, subsample=1,
                                                                   max_depth=3, verbose=2))
    ])

    print('Pipeline created')
    if cross_val:
        print('Five-fold cross-validation')
        five_fold_cross_val(pipeline, articles, labels)
    print('Fitting classifier')
    pipeline.fit(articles, labels)

    return pipeline