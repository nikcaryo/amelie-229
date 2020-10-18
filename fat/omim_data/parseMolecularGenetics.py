#! /usr/bin/env python

import lxml.etree as ET
import sys
import re

if __name__ == "__main__":
    omim_id = sys.argv[1].split('/')[-1].split('.')[0]
    with open(sys.argv[1]) as f:
        doc = ET.parse(f)
    mgs = doc.xpath(".//textSection[contains(textSectionName, 'molecularGenetics')]/textSectionContent")
    pmid_pattern = re.compile('\{([0-9]+):.*?\}')
    references = []
    for mg in mgs:
        content = ''.join(mg.itertext())
        # replaced = re.findall('((previously|originally) (described|reported|studied) by [^{]*{[^}]*}[^,\.]*[\.,])', content)
        # for r in replaced:
        #     print(r[0])
        content = re.sub('((previously|originally) (described|reported|studied) by [^{]*{[^}]*}[^,\.]*[\.,])', '', content)
        references.extend([f for f in pmid_pattern.findall(content)])
    references = set(references)
    pmids = set()
    refs = doc.xpath('//reference')
    for ref in refs:
        pubmedIDs = ref.xpath("pubmedID")
        if len(pubmedIDs) == 0:
            continue
        assert len(pubmedIDs) == 1
        referenceNumber = ''.join(ref.xpath("referenceNumber")[0].itertext())
        if referenceNumber in references:
            pmids.add(''.join(pubmedIDs[0].itertext()))
    for pmid in pmids:
        print("%s\t%s" % (omim_id, pmid))
