#! /usr/bin/env python3

import libxml2
import sys
import re

if __name__ == "__main__":
  doc = libxml2.parseFile(sys.argv[1])
  ctxt = doc.xpathNewContext()
  refs = ctxt.xpathEval("//reference")
  for ref in refs:
    ctxt.setContextNode(ref)
    pubmedIDs = ctxt.xpathEval("pubmedID")
    if len(pubmedIDs) == 0:
      continue
    assert len(pubmedIDs) == 1
    omim_id = sys.argv[1].split('/')[1].split('.')[0]
    print(omim_id + "\t" + pubmedIDs[0].getContent())
