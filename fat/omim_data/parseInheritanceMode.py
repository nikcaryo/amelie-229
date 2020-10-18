#! /usr/bin/env python3

from lxml import etree
import sys
import re

if __name__ == "__main__":
    for line in sys.stdin:
        omim_id = line.strip().split('/')[-1].split('.')[0]
        doc = etree.parse(open(line.strip()))
        mgs = doc.xpath("//phenotypeInheritance")
        texts = []
        for mg in mgs:
            texts.append(str(''.join([x for x in mg.itertext()])))
        inheritance_modes = set(['DOMINANT' if 'dominant' in x else ('RECESSIVE' if 'recessive' in x else 'UNKNOWN') for x in texts])
        if len(inheritance_modes) > 1:
            inheritance_mode = 'UNKNOWN'
        elif len(inheritance_modes) == 0:
            inheritance_mode = 'UNKNOWN'
        else:
            inheritance_mode = list(inheritance_modes)[0]
        print("%s\t%s" % (omim_id, inheritance_mode))
