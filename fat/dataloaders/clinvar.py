
import csv
import io

import pkg_resources

from collections import defaultdict

from amelie.tuples import VCFCoords


def parent_package_name():
    return '.'.join(__name__.split('.')[:-2])


class Clinvar:
    def __init__(self):
        self.mapping = defaultdict(set)
        with pkg_resources.resource_stream(parent_package_name(),
                                           'avada_data/'
                                           'clinvar20160705_pathogenic_variants_realigned_with_pmids.tsv') \
                as clinvar_data:
            with io.TextIOWrapper(clinvar_data) as clinvar_data_text:
                for row in csv.reader(clinvar_data_text, delimiter='\t'):
                    self.mapping[int(row[4])].add(VCFCoords(chrom=row[0],
                                                            pos=int(row[1]),
                                                            ref=row[2],
                                                            alt=row[3]))

    def get_labeled_pmids(self):
        return self.mapping.keys()

    def get_dna_changes(self, pmid):
        return self.mapping.get(pmid, set())