#! /usr/bin/env python3

import sys

def read_ensembl_genes(filename):
    rv = []
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == 0 and line.startswith('Ensembl'):
                continue
            line = [str(l) for l in line.strip().split('\t')]
            rv.append(line)
    return rv

def read_hgnc_symbols(filename):
    rv = {}
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == 0 and line.startswith('HGNC'):
                continue
            line = line.strip().split('\t')
            while len(line) < 4:
                line.append('')
            line[3] = [s.strip() for s in line[3].split(', ') if s.strip() != '']
            rv[line[0]] = tuple(line[1:])
    return rv

if __name__ == "__main__":
    ensembl_genes = read_ensembl_genes(sys.argv[1])
    hgnc_symbols = read_hgnc_symbols(sys.argv[2])
    for ensembl_gene in ensembl_genes:
        if len(ensembl_gene) >= 2 and ensembl_gene[2] in hgnc_symbols:
            hgnc_entry = hgnc_symbols[ensembl_gene[2]]
            hgnc_name = hgnc_entry[0]
            hgnc_synonyms = hgnc_entry[2]
            if hgnc_name != ensembl_gene[1]:
                sys.stderr.write("Warning: ensembl gene id %s: ensembl gene name %s != hgnc symbol %s\n" % (ensembl_gene[0], ensembl_gene[1], hgnc_name))
            print("%s\t%s\t%s\t%s" % (ensembl_gene[0], hgnc_name, hgnc_name, 'HGNC_SYMBOL'))
            for synonym in hgnc_synonyms:
                print("%s\t%s\t%s\t%s" % (ensembl_gene[0], hgnc_name, synonym, 'HGNC_SYNONYM'))
        else:
            sys.stderr.write("Skipping ensembl gene id %s: doesn't have HGNC symbol\n" % ensembl_gene[0])
