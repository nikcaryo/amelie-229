from collections import defaultdict
import urllib.request
import os.path
import gzip
from .. import tuples
from . import gene_data
import constdb
import pickle
import io
import csv

pubtator_path = 'ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator/gene2pubtator.gz'

#     pubtator.download_pubtator(args.out_dir + '/pubtator.db')


def download_pubtator(output_path):
    ncbi2eid = gene_data.load_entrez_to_ensembl()
    eid2symbol = gene_data.ensembl_id_map_to_gene_symbol()

    pmid2genes = defaultdict(lambda: defaultdict(set))

    with urllib.request.urlopen(pubtator_path) as pubtator_file:
        with gzip.open(pubtator_file) as pubtator_bytes:
            with io.TextIOWrapper(pubtator_bytes) as pubtator_text:
                for row in csv.DictReader(pubtator_text, delimiter='\t'):
                    pmid = int(row['PMID'])
                    mention_list = row['Mentions'].split("|")

                    if row['NCBI_Gene'] in ncbi2eid:
                        eid = ncbi2eid[row['NCBI_Gene']]
                    else:
                        continue
                    if eid not in eid2symbol:
                        continue

                    for mention in mention_list:
                        mention = gene_data.normalize_gene_name(mention.strip().split(), 'PROTEIN_GENE2PUBTATOR')
                        pmid2genes[pmid][mention].add(tuples.GeneName(eid, eid2symbol[eid], 'PROTEIN_GENE2PUBTATOR'))

    with constdb.create(output_path) as db:
        for pmid, mention_dict in pmid2genes.items():
            mention_pickle = pickle.dumps(mention_dict)
            db.add(pmid, mention_pickle)


class Pubtator:
    def __init__(self, path):
        self.db = constdb.read(path, mmap=False)

    def get_pmid2gene(self, pmid):
        pmid_data = self.db.get(pmid)

        if pmid_data is None:
            return {}
        else:
            return pickle.loads(pmid_data)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.db.close()
        return False