import csv
import pkg_resources

import io
from .clinvar import parent_package_name

def get_gwas_pmids():
    result = set()
    with pkg_resources.resource_stream(parent_package_name(), 'gwas_data/gwas_catalog_v1.0-studies_r2018-04-10.tsv') as gwas_data:
        with io.TextIOWrapper(gwas_data) as gwas_file:
            for row in csv.DictReader(gwas_file, delimiter='\t'):
                result.add(int(row['PUBMEDID']))

    return result