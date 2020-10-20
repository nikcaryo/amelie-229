#! /bin/zsh

wget 'http://purl.obolibrary.org/obo/hp.obo' -O hpo_ontology.obo

# MeSH (downloaded October 15, 2018):

wget 'ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/asciimesh/c2018.bin' 
wget 'ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/asciimesh/d2018.bin'

# SNOMED was downloaded from: https://download.nlm.nih.gov/mlb/utsauth/USExt/SnomedCT_USEditionRF2_PRODUCTION_20180901T120000Z.zip on October 15, 2018
#
# UMLS 2018AA downloaded from https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html on October 15, 2018

# stats about xrefs in HPO:

cat hpo_ontology.obo | egrep '^xref' | cut -d ':' -f 2 | sort | uniq -c | sort -k1,1n
cat hpo_ontology.obo | egrep '^xref' | cut -d ':' -f 2 | sort | uniq -c | awk '{print $1}' | datamash sum 1

# UMLS:

zcat umls/2018AA-full/2018AA/META/MRCONSO.RRF.a*.gz | cut -d '|' -f 1,15 | tr '|' '\t' |
  sort -u |
  sort -t $'\t' -k1,1 > umls_terms.tsv

cat hpo_ontology.obo | grep xref | grep UMLS | cut -d ':' -f 3 | sort -u > hpo_umls.txt

join -t $'\t' -1 1 -2 1 umls_terms.tsv hpo_umls.txt > hpo_umls_terms.tsv

# SNOMED:

cut -f 5,8 snomed/SnomedCT_USEditionRF2_PRODUCTION_20180901T120000Z/Full/Terminology/sct2_Description_Full-en_US1000124_20180901.txt | 
  sort -u | sort -t $'\t' -k1,1 > snomed_terms.tsv

cat hpo_ontology.obo | grep xref | grep SNOMED | cut -d ':' -f 3 | sort -u > hpo_snomed.txt

join -t $'\t' -1 1 -2 1 snomed_terms.tsv hpo_snomed.txt > hpo_snomed_terms.tsv

# MeSH:

./mesh/parse_mesh_c.py mesh/c2018.bin > mesh_terms.tsv
./mesh/parse_mesh_d.py mesh/d2018.bin >> mesh_terms.tsv

cat hpo_ontology.obo | grep xref | grep MSH | cut -d ':' -f 3 | sort -u > hpo_mesh.txt

join -t $'\t' -1 1 -2 1 mesh_terms.tsv hpo_mesh.txt > hpo_mesh_terms.tsv
