#! /usr/bin/env python3

import xml.sax
import sys
from collections import defaultdict

class UniprotHandler(xml.sax.ContentHandler):
    
    def __init__(self):
        self.cur_content = ""
        self.mode = ""
        self.cur_accession_ids = []
        self.cur_names = {}

    def characters(self, content):
        self.cur_content += content

    def startElement(self, name, attrs):
        sys.stdout.flush()
        if name == "accession" or \
            name == "fullName" or \
            name == "shortName" or \
            name == "name":
            self.cur_content = ""
        if name == "protein":
            self.mode = "protein"
        if name == "gene":
            self.mode = "gene"

    def endElement(self, name):
        if name == "accession":
            self.cur_accession_ids.append(self.cur_content)
        if name == "protein":
            self.mode = ""
        if name == "gene":
            self.mode = ""
        if name == "fullName" and self.mode == "protein":
            self.cur_names[self.cur_content] = "PROTEIN_FULL_NAME"
        if name == "shortName" and self.mode == "protein":
            self.cur_names[self.cur_content] = "PROTEIN_SHORT_NAME"
        if name == "name" and self.mode == "protein":
            self.cur_names[self.cur_content] = "PROTEIN_NAME"
        if name == "fullName" and self.mode == "gene":
            self.cur_names[self.cur_content] = "UNIPROT_GENE_NAME"
        if name == "shortName" and self.mode == "gene":
            self.cur_names[self.cur_content] = "UNIPROT_GENE_NAME"
        if name == "name" and self.mode == "gene":
            self.cur_names[self.cur_content] = "UNIPROT_GENE_NAME"
        if name == "entry":
            for accession_id in self.cur_accession_ids:
                for name in self.cur_names:
                    print("%s\t%s\t%s" % (accession_id, name, self.cur_names[name]))
                    sys.stdout.flush()
            self.cur_accession_ids = []
            self.cur_names = {}

if __name__ == "__main__":
    parser = xml.sax.make_parser()
    u = UniprotHandler()
    parser.setContentHandler(u)
    with open(sys.argv[1]) as f:
        parser.parse(f)
