
from collections import defaultdict
import pkg_resources


def load_ensembl_to_refseq():
    refseq_to_ensembl_tsv = pkg_resources.resource_string(__name__, 'gene_data/nm_to_ensembl_to_np.tsv').decode('utf8')
    rv = defaultdict(lambda: set())
    for line in refseq_to_ensembl_tsv.split('\n'):
        if len(line) == 0:
            continue
        nm_id, ensembl_id, np_id = line.split('\t')
        if nm_id != "":
            rv[ensembl_id].add(nm_id)
        if np_id != "":
            rv[ensembl_id].add(np_id)
    return rv


class RefseqExtractor:
    def __init__(self):
        self.ensembl_to_refseq = load_ensembl_to_refseq()

    def extract(self, words):
        rv = defaultdict(lambda: set())
        all_words = words.split()
        for ensembl_id, refseq_ids in self.ensembl_to_refseq.items():
            for refseq_id in refseq_ids:
                if refseq_id in all_words:
                    rv[ensembl_id].add(refseq_id)
        return rv
