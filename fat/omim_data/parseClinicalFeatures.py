#! /usr/bin/env python3

from lxml import etree
import sys
import re
from nltk import word_tokenize

sys.path.append('../../src')

from text import util
from text.extractors.pheno_extractor import extract_candidate_mentions, load_pheno_terms
from text import extractor_util

PHENOS, PHENO_SETS = load_pheno_terms()
CANON_PHENO_NAMES = util.load_canon_pheno_names()

if __name__ == "__main__":
  for line in sys.stdin:
    omim_id = line.strip().split('/')[-1].split('.')[0]
    doc = etree.parse(open(line.strip()))
    mgs = doc.xpath("//textSection[contains(textSectionName, 'clinicalFeatures')]/textSectionContent")
    text = ""
    for mg in mgs:
      text += str("".join([x for x in mg.itertext()]))
    paragraphs = text.split('\n\n')
    for i, paragraph in enumerate(paragraphs):
      words = word_tokenize(paragraph)
      mentions = extract_candidate_mentions(words, PHENOS, PHENO_SETS, CANON_PHENO_NAMES)
      print("\n".join(omim_id + "\t" + str(i) + "\t" + str(x.wordidxs[0]) + "\t" + x.hpo_id for x in sorted(mentions)))
      sys.stdout.flush()
